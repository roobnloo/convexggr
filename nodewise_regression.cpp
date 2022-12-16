#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

VectorXd softThreshold(const VectorXd &x, double lambda)
{
    MatrixXd thresh = MatrixXd::Zero(x.rows(), 2);
    thresh.col(0) = x.array().abs() - lambda;
    return x.array().sign() * thresh.array().rowwise().maxCoeff();
}

double softThreshold(double x, double lambda)
{
    int signx = (double(0) < x) - (x < double(0));
    return signx * std::max(abs(x) - lambda, 0.0);
}

// [[Rcpp::export]]
VectorXd applyRidgeUpdate(
    const VectorXd &v, VectorXd &residual, const MatrixXd &U, double regmean)
{
    int p = U.rows();
    auto regdiag = VectorXd::Constant(p, regmean).asDiagonal();
    residual += U * v;
    MatrixXd UtUreg = U.transpose() * U;
    UtUreg += regdiag;
    VectorXd v_update = UtUreg.llt().solve(U.transpose() * residual);
    residual -= U * v_update;

    return v_update;
}

// [[Rcpp::export]]
VectorXd applyL1Update(
    const VectorXd &b, VectorXd &residual, const MatrixXd &X, double penalty)
{
    int n = X.rows();
    int p = b.rows();

    if (residual.rows() != n || X.cols() != p)
    {
        stop("Dimension mismatch during L1 update step!");
    }

    VectorXd v_update(p);
    for (int i = 0; i < b.rows(); ++i)
    {
        auto uslice = X.col(i);
        VectorXd presid(residual);
        presid += b(i) * uslice;
        v_update(i) = softThreshold(b(i) + residual.dot(uslice) / n, penalty);
        residual = presid - v_update(i) * uslice;
    }

    return v_update;
}

double objective(
    const VectorXd &residual, const VectorXd &mean_coef, const MatrixXd &beta,
    double regmean, double lambda, double asparse)
{
    double quad_loss = residual.squaredNorm() / (2 * residual.rows());
    double mean_obj = mean_coef.squaredNorm();
    double lasso_obj = beta.col(0).lpNorm<1>();
    double group_lasso_obj = 0;
    for (int i = 1; i < beta.cols(); ++i)
    {
        VectorXd bcol(beta.col(i));
        lasso_obj += bcol.lpNorm<1>();
        group_lasso_obj += bcol.norm();
    }
    return quad_loss +
           regmean * mean_obj +
           asparse * lambda * lasso_obj +
           (1 - asparse) * lambda * group_lasso_obj;
}

double objective_sgl(
    const VectorXd &residual, const VectorXd &v, double lambda, double asparse)
{
    int n = residual.rows();
    return residual.squaredNorm() / (2 * n) +
           asparse * lambda * v.lpNorm<1>() +
           (1 - asparse) * lambda * v.norm();
}

VectorXd sglUpdateStep(
    const VectorXd &center, const VectorXd &grad,
    double step, double lambda, double asparse)
{
    VectorXd thresh = softThreshold(center - step * grad,
                                    step * asparse * lambda);
    double threshnorm = thresh.norm();
    if (threshnorm == 0)
    {
        return VectorXd::Zero(center.rows());
    }
    double maxterm = std::max(0.0, 1 - step * (1 - asparse) * lambda / threshnorm);
    return maxterm * thresh;
}

// [[Rcpp::export]]
VectorXd applySparseGLUpdate(
    const VectorXd &beta_grp, VectorXd &residual,
    const VectorXd &covariate, const MatrixXd &response,
    double lambda, double asparse, int maxit, double tol)
{
    int n = response.rows();
    MatrixXd intx = response.array().colwise() * covariate.array();
    residual += intx * beta_grp;
    VectorXd threshold = softThreshold(
        intx.transpose() * residual / n, asparse * lambda);

    // If this subgradient condition holds, the entire group should be zero
    if (threshold.norm() <= (1 - asparse) * lambda)
    {
        return VectorXd::Zero(beta_grp.rows());
    }

    double step_size = 1;
    VectorXd beta_update(beta_grp);
    VectorXd center_new(beta_grp);
    double obj = INFINITY;
    double objnew;
    for (int i = 0; i < maxit; ++i)
    {
        objnew = objective_sgl(
            residual - intx * beta_update, beta_update, lambda, asparse);
        if (abs(objnew - obj) < tol)
        {
            break;
        }
        if (i == maxit)
        {
            warning(
                "Inner loop within SGL descent exceeded max %i iterations",
                maxit);
        }
        obj = objnew;
        VectorXd center_old = center_new;
        VectorXd grp_fit = intx * beta_update;
        VectorXd grad = -1 * intx.transpose() * (residual - grp_fit) / n;

        // Optimize the step size
        double lhs = INFINITY;
        double rhs = 0;
        double quad_loss_old = (residual - grp_fit).squaredNorm() / (2 * n);
        while (true)
        {
            center_new = sglUpdateStep(
                beta_update, grad, step_size, lambda, asparse);
            VectorXd centerdiff = center_new - beta_update;
            rhs = quad_loss_old + grad.dot(centerdiff) +
                  centerdiff.squaredNorm() / (2 * step_size);
            lhs = (residual - intx * center_new).squaredNorm() / (2 * n);
            if (lhs <= rhs)
                break;
            step_size *= 0.8;
        }

        // Nesterov momentum step
        beta_update = center_old +
                      ((i + 1.0) / (i + 4)) * (center_new - center_old);
    }

    residual -= intx * beta_update;
    return beta_update;
}

// TODO
double estimateVariance(
    const VectorXd &residual, const VectorXd &mean_coef, const MatrixXd &beta)
{
    return 1;
}

// [[Rcpp::export]]
List nodewiseRegression(
    VectorXd y, MatrixXd response, MatrixXd covariates,
    double lambda, double asparse, double regmean,
    VectorXd initbeta, VectorXd initgamma, int maxit = 1000, double tol = 1e-5)
{
    int p = response.cols() + 1;
    int q = covariates.cols();
    y = y.array() - y.mean();

    // MatrixXd beta(initbeta);
    // beta.conservativeResize(p-1, q+1);
    // VectorXd gamma(initgamma);
    // TODO: fix initialization
    MatrixXd beta = MatrixXd::Zero(p - 1, q + 1);
    VectorXd gamma = VectorXd::Zero(q);

    VectorXd residual = y - covariates * gamma - response * beta.col(0);
    for (int i = 0; i < q; ++i)
    {
        MatrixXd intxgrp =
            response.array().colwise() * covariates.col(i).array();
        residual -= intxgrp * beta.col(i + 1);
    }

    NumericVector objval(maxit + 1);
    objval[0] = objective(residual, gamma, beta, regmean, lambda, asparse);
    Rcout << "First objective: " << objval[0] << std::endl;

    for (int i = 0; i < maxit + 1; ++i)
    {
        gamma = applyRidgeUpdate(gamma, residual, covariates, regmean);

        beta.col(0) = applyL1Update(
            beta.col(0), residual, response, lambda * asparse);

        for (int j = 0; j < q; ++j)
        {
            beta.col(j + 1) = applySparseGLUpdate(
                beta.col(j + 1), residual, covariates.col(j), response,
                lambda, asparse, maxit, tol);
        }

        objval[i + 1] = objective(residual, gamma, beta, regmean, lambda, asparse);

        if (i + 1 == 1)
        {
            Rcout << "Second objective: " << objval[i + 1] << std::endl;
        }

        if (abs(objval[i + 1] - objval[i]) < tol)
        {
            Rcout << "Breaking after (" << (i + 1) << ") iterations.";
            objval = objval[Rcpp::Range(0, i + 1)];
            break;
        }
    }
    if (objval.length() >= maxit)
    {
        warning("Maximum iterations exceeded!");
    }

    VectorXd betahat = VectorXd::Zero((p - 1) * (q + 1));
    for (int i = 0; i < q + 1; ++i)
    {
        betahat.segment(i, p - 1) = beta.col(i);
    }

    return List::create(
        Named("beta") = betahat,
        Named("gamma") = gamma,
        Named("objval") = objval,
        Named("resid") = residual);
}

// TODO: handle exceptions