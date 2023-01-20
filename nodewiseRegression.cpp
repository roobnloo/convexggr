#include <Rcpp.h>
#include <RcppEigen.h>
#include <limits>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// For easy debugging
/*
void print(VectorXd v)
{
    std::cout << v << std::endl;
}

void print(MatrixXd m, int c)
{
    std::cout << m.col(c) << std::endl;
}
*/

double softThreshold(double x, double lambda)
{
    double diff = std::abs(x) - lambda;
    if (diff < 0)
    {
        return 0;
    }
    int signx;
    if (std::abs(x) <= std::numeric_limits<double>::epsilon())
    {
        signx = 0;
    }
    else
    {
        signx = (double(0) < x) - (x < double(0));
    }
    return signx * diff;
}

VectorXd softThreshold(const VectorXd &x, double lambda)
{
    VectorXd thresh(x.rows());
    for (int i = 0; i < x.rows(); ++i)
    {
        thresh(i) = softThreshold(x(i), lambda);
    }
    return thresh;
}

VectorXd applyRidgeUpdate(
    const VectorXd &v, VectorXd &residual, const MatrixXd &U, double regmean)
{
    if (U.cols() != v.rows())
    {
        stop("Dimension mismatch during ridge update step!");
    }

    int p = U.rows();
    auto regdiag = VectorXd::Constant(p, regmean).asDiagonal();
    residual += U * v;
    MatrixXd UtUreg = U.transpose() * U;
    UtUreg += regdiag;
    VectorXd v_update = UtUreg.llt().solve(U.transpose() * residual);
    residual -= U * v_update;

    return v_update;
}

VectorXd applyL1Update(
    const VectorXd &b, VectorXd &residual, const MatrixXd &X, double penalty)
{
    int n = X.rows();
    int p = b.rows();

    if (residual.rows() != n || X.cols() != p)
    {
        stop("Dimension mismatch during L1 update step!");
    }

    VectorXd bnext(b);
    for (int i = 0; i < b.rows(); ++i)
    {
        auto xslice = X.col(i);
        residual += bnext(i) * xslice;
        double xscale = xslice.squaredNorm() / n;
        bnext(i) = softThreshold(residual.dot(xslice) / n, penalty) / xscale;
        residual -= bnext(i) * xslice;
    }

    return bnext;
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
        lasso_obj += beta.col(i).lpNorm<1>();
        group_lasso_obj += beta.col(i).norm();
    }
    double result = quad_loss +
                    regmean * mean_obj +
                    asparse * lambda * lasso_obj +
                    (1 - asparse) * lambda * group_lasso_obj;
    return result;
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
    VectorXd thresh = softThreshold(
        center - step * grad, step * asparse * lambda);
    double threshnorm = thresh.norm();
    if (threshnorm <= step * (1 - asparse) * lambda)
    {
        return VectorXd::Zero(center.rows());
    }
    double normterm = 1 - step * (1 - asparse) * lambda / threshnorm;
    if (normterm < 0)
    {
        return VectorXd::Zero(center.rows());
    }
    return normterm * thresh;
}

VectorXd applySparseGLUpdate(
    const VectorXd &beta_grp, VectorXd &residual, const MatrixXd &intx,
    double lambda, double asparse, int maxit, double tol)
{
    int n = intx.rows();
    if (residual.rows() != n || intx.cols() != beta_grp.rows())
    {
        stop("Dimension mismatch during sparse group lasso update step!");
    }

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
        if (std::abs(objnew - obj) < tol)
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

// TODO: is the variance correct?
double estimateVariance(
    const VectorXd &residual, const VectorXd &gamma, const MatrixXd &beta)
{
    int numNonZero = (beta.array().abs() > 0).count();
    return residual.squaredNorm() / (residual.rows() - numNonZero);
}

NumericVector getLambda(
    NumericVector inlambda, int nlambda, double lambdaFactor,
    const VectorXd &y, const std::vector<MatrixXd> &intxs)
{
    if (inlambda.size() > 0)
    {
        NumericVector inSorted = inlambda.sort(true);
        if (inSorted[inSorted.size() - 1] <= 0)
        {
            stop("Provided lambda terms are not all positive!");
        }
        return inSorted;
    }
    double lamdaMax = 0;
    for (MatrixXd intx : intxs)
    {
        double mc = (intx.transpose() * y).array().abs().maxCoeff();
        if (mc > lamdaMax)
        {
            lamdaMax = mc;
        }
    }

    NumericVector loglinInterp(nlambda);
    double delta = log(lambdaFactor) / (nlambda - 1);
    for (int i = 0; i < nlambda; ++i)
    {
        loglinInterp[i] = lamdaMax * exp(i * delta);
    }
    return loglinInterp;
}

struct RegressionResult {
    VectorXd resid;
    double varhat;
    double objval;
};

RegressionResult nodewiseRegressionInit(
    const VectorXd &y, const MatrixXd &response, const MatrixXd &covariates,
    const std::vector<MatrixXd> &intxs,
    VectorXd &gamma, MatrixXd &beta, // initial guess
    double lambda, double asparse, double regmean,
    int maxit, double tol, bool verbose)
{
    int p = response.cols() + 1;
    int q = covariates.cols();
    beta.resize(p-1, q+1);

    VectorXd residual = y - covariates * gamma - response * beta.col(0);
    for (int i = 0; i < (int) intxs.size(); ++i)
    {
        residual -= intxs[i] * beta.col(i + 1);
    }

    NumericVector objval(maxit + 1);
    objval[0] = objective(residual, gamma, beta, regmean, lambda, asparse);

    for (int i = 0; i < maxit; ++i)
    {
        gamma = applyRidgeUpdate(gamma, residual, covariates, regmean);

        beta.col(0) = applyL1Update(
            beta.col(0), residual, response, lambda * asparse);

        for (int j = 0; j < q; ++j)
        {
            beta.col(j + 1) = applySparseGLUpdate(
                beta.col(j + 1), residual, intxs[j],
                lambda, asparse, maxit, tol);
        }

        objval[i + 1] = objective(residual, gamma, beta, regmean, lambda, asparse);
        // if (verbose)
        //     std::cout << "Iteration: " << i << ":: obj:" << objval[i+1] << std::endl;
        
        if (i > 4 && std::abs(objval[i + 1] - objval[i - 1]) < 1e-20
            && std::abs(objval[i] - objval[i - 2]) < 1e-20)
        {
            // std::cout << "potential oscillation" << std::endl;
            stop("Potential oscillation!");
        }

        if (std::abs(objval[i + 1] - objval[i]) < tol)
        {
            objval = objval[Rcpp::Range(0, i + 1)];
            break;
        }
    }
    if (verbose)
        std::cout << "Finished in " << objval.length() << " iterations" << std::endl;
    if (objval.length() == maxit + 1)
    {
        std::cout << "Maximum iterations exceeded!" << std::endl;
        warning("Maximum iterations exceeded!");
    }

    beta.resize((p-1)*(q+1), 1);
    double varhat = estimateVariance(residual, gamma, beta);

    return RegressionResult {
        residual, varhat, objval[objval.size() - 1]
    };
}

// [[Rcpp::export]]
List nodewiseRegression(
    VectorXd y, MatrixXd response, MatrixXd covariates,
    double asparse, double regmean,
    NumericVector lambdas = NumericVector::create(),
    int nlambda = 100, double lambdaFactor = 1e-4,
    int maxit = 1000, double tol = 1e-8, bool verbose = false)
{
    int p = response.cols() + 1;
    int q = covariates.cols();

    if (response.rows() != covariates.rows() || y.rows() != response.rows())
    {
        stop("Responses and covariates must have the same number of observations!");
    }
    if (asparse < 0 || asparse > 1 || regmean < 0)
    {
        stop("Penalty terms are out of range!");
    }
    if (maxit <= 0 || tol <= 0)
    {
        stop("Maximium iteration and/or numerical tolerance are out of range!");
    }

    std::vector<MatrixXd> intxs(q);
    for (int i = 0; i < q; ++i)
    {
        MatrixXd intx =
            response.array().colwise() * covariates.col(i).array();
        intxs[i] = intx;
    }

    lambdas = getLambda(lambdas, nlambda, lambdaFactor, y, intxs);
    nlambda = lambdas.size(); // a bit of a hack for user-provided lambdas

    MatrixXd gammaFull(q, nlambda);
    MatrixXd betaFull((p-1)*(q+1), nlambda);
    VectorXd varhatFull(nlambda);
    MatrixXd residualFull(y.rows(), nlambda);
    VectorXd objectiveFull(nlambda);

    RegressionResult regResult;
    MatrixXd beta(MatrixXd::Zero((p-1)*(q+1), 1));
    VectorXd gamma(VectorXd::Zero(q));
    for (int i = 0; i < nlambda; ++i)
    {
        if (verbose)
            std::cout << "Regression with lambda index " << i << std::endl;
        regResult = nodewiseRegressionInit(
            y, response, covariates, intxs, gamma, beta,
            lambdas[i], asparse, regmean, maxit, tol, verbose);

        // Use gamma and beta as initializers for next lambda (warm-starts)
        gammaFull.col(i) = gamma;
        betaFull.col(i) = beta;
        residualFull.col(i) = regResult.resid;
        varhatFull(i) = regResult.varhat;
        objectiveFull(i) = regResult.objval;
    }

    return List::create(
        Named("beta") = betaFull,
        Named("gamma") = gammaFull,
        Named("varhat") = varhatFull,
        Named("objval") = objectiveFull,
        Named("resid") = residualFull,
        Named("lambdas") = lambdas);
}