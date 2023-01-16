library(abind)
library(Rcpp)
source("utils.R")
source("cv_cggr_node.R")
# source("node_strategy_sparsegl.R")
sourceCpp("nodewiseRegression.cpp")

#' @param responses n x d matrix of responses
#' @param U n x p matrix of covariates
#' @param regmean small ridge penalty on mean vector
#' @param asparse SGL mixture penalty
cggr <- function(responses, covariates, asparse,
                 regmean = 0.01, nlambda = 100, lambdafactor = 1e-4,
                 maxit = 3e+06, tol = 1e-8, nfolds = 5, verbose = FALSE) {
  stopifnot(is.matrix(responses), is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            regmean > 0)

  p <- ncol(responses)
  q <- ncol(covariates)
  n <- nrow(responses)

  # Initial run to get lambdas for each response
  lambdas <- matrix(nrow = nlambda, ncol = p)
  beta <- array(dim = c((p - 1) * (q + 1), nlambda, p))
  gamma <- array(dim = c(q, nlambda, p))
  varhat <- matrix(nrow = nlambda, ncol = p)
  resid <- array(dim = c(n, nlambda, p))
  objval <- matrix(nrow = nlambda, ncol = p)

  for (node in seq_len(p)) {
    nodereg <- nodewiseRegression(
      responses[, node], responses[, -node], covariates, asparse, regmean,
      nlambda = nlambda, lambdaFactor = lambdafactor, maxit = maxit, tol = tol)

    lambdas[, node] <- nodereg["lambdas"][[1]]
    beta[, , node] <- nodereg["beta"][[1]]
    gamma[, , node] <- nodereg["gamma"][[1]]
    varhat[, node] <- nodereg["varhat"][[1]]
    resid[, , node] <- nodereg["resid"][[1]]
    objval[, node] <- nodereg["objval"][[1]]

    if (verbose)
      print(paste("Finished initial run for node", node))
  }

  if (verbose)
    print("Finished initial run")

  cv_lambda_idx <- vector(length = p)
  cv_mse <- matrix(nrow = p, ncol = nlambda)
  for (node in seq_len(p)) {
    if (verbose)
      print(paste("CV for node", node))
    cv_result <- cv_cggr_node(
      node, responses, covariates,
      lambdas[, node], asparse, regmean, maxit, tol, nfolds, verbose)
    cv_mse[node, ] <- cv_result$cv_mse
    cv_lambda_idx[node] <- which.min(cv_result$cv_mse)
  }

  if (verbose)
    print("Finished cross validating all nodes")

  return(list(lambda = lambdas,
              beta = beta,
              objval = objval,
              cv_mse = cv_mse,
              cv_lambda = cv_lambda_idx))

  # # Initialize mean matrix
  # ghat_mx <- matrix(nrow = p, ncol = q)

  # # Initialize covariate array
  # # Includes the population matrix, hence +1
  # bhat_tens <-  array(0, dim = c(p, p, q + 1))

  # # Estimated variances
  # varhat <- vector(length = p)

  # for (node in seq_len(p)) {
  #   result <- nodereg_strat(
  #     node, X, U,
  #     lambda_min[node], asparse, regmean,
  #     initbeta = rep(0, (p - 1) * (q + 1)), initgamma = rep(0, q))

  #   # Do we need a scaling factor here? I think so.
  #   varhat[node] <- result$varhat
  #   bhat_tens[node, -node, ] <- -result$beta / varhat[node]
  #   ghat_mx[node, ] <- result$gamma
  #   if (verbose) {
  #     print(paste("Finished regression for node", node,
  #                 "with obj", result$objval[length(result$objval)]))
  #   }
  # }

  # bhat_symm <- abind(apply(bhat_tens, 3, symmetrize, simplify = F), along = 3)

  # # Returns the estimated precision matrix of the ith observation
  # precision <- function(i) {
  #   Theta <- apply(bhat_symm, c(1, 2), \(b) b %*% c(1, U[i, ]))
  #   diag(Theta) <- 1/varhat
  #   # omega <- -1 * diag(1/varhat) %*% Theta
  #   return(Theta)
  # }

  # # Returns the estimated mean vector of the ith observation
  # mean <- function(i) {
  #   prec <- precision(i)
  #   mu <- solve(prec, diag(1 / varhat) %*% ghat_mx %*% U[i, ])
  #   return(mu)
  # }

  # return(list(ghat = ghat_mx,
  #             bhat = bhat_symm,
  #             bhat_asym = bhat_tens,
  #             varhat = varhat,
  #             precision = precision,
  #             mean = mean))
}
