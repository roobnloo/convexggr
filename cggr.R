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
                 maxit = 1000, tol = 1e-8, nfolds = 5,
                 verbose = FALSE, parallel = FALSE) {

  stopifnot(is.matrix(responses), is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            regmean > 0)

  p <- ncol(responses)
  q <- ncol(covariates)
  n <- nrow(responses)
  bveclength <- (p - 1) * (q + 1)

  # Initial run to get lambdas for each response
  lambdas <- matrix(nrow = nlambda, ncol = p)
  beta <- array(dim = c(bveclength, nlambda, p))
  gamma <- array(dim = c(q, nlambda, p))
  varhat <- matrix(nrow = nlambda, ncol = p)
  resid <- array(dim = c(n, nlambda, p))
  objval <- matrix(nrow = nlambda, ncol = p)

  nodewise <- function(node) {
    nodereg <- nodewiseRegression(
      responses[, node], responses[, -node], covariates, asparse, regmean,
      nlambda = nlambda, lambdaFactor = lambdafactor, maxit = maxit, tol = tol)
    if (verbose)
      print(paste("Finished initial run for node", node))
    return(list(
      lambdas = nodereg["lambdas"][[1]],
      beta = nodereg["beta"][[1]],
      gamma = nodereg["gamma"][[1]],
      varhat = nodereg["varhat"][[1]],
      resid = nodereg["resid"][[1]],
      objval = nodereg["objval"][[1]]
    ))
  }

  if (parallel) {
    print("Running in parallel...")
    tictoc::tic()
    reg_result <- parallel::mclapply(seq_len(p), nodewise)
    tictoc::toc()
  } else {
    tictoc::tic()
    reg_result <- lapply(seq_len(p), nodewise)
    tictoc::toc()
  }
  for (node in seq_len(p)) {
    lambdas[, node] <- reg_result[[node]]$lambdas
    beta[, , node] <- reg_result[[node]]$beta
    gamma[, , node] <- reg_result[[node]]$gamma
    varhat[, node] <- reg_result[[node]]$varhat
    resid[, , node] <- reg_result[[node]]$resid
    objval[, node] <- reg_result[[node]]$objval
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

  # return(list(lambda = lambdas,
  #             beta = beta,
  #             objval = objval,
  #             cv_mse = cv_mse,
  #             cv_lambda = cv_lambda_idx))

  ghat_select <- cbind(rep(seq(q), times = p),
                       rep(cv_lambda_idx, each = q),
                       rep(seq(p), each = q))
  ghat_mx <- t(matrix(gamma[ghat_select], nrow = q, ncol = p))

  varhat <- varhat[cbind(cv_lambda_idx, seq(p))]

  # Includes the population matrix, hence +1
  bhat_select <- cbind(rep(seq(bveclength), times = p),
                       rep(cv_lambda_idx, each = bveclength),
                       rep(seq(p), each = bveclength))
  bhat_mx <- matrix(beta[bhat_select], nrow = bveclength, ncol = p)
  bhat_tens <-  array(0, dim = c(p, p, q + 1))
  for (i in seq_len(p)) {
    bhat_tens[i, -i, ] <- -bhat_mx[, i] / varhat[i]
  }

  bhat_symm <- abind(
    apply(bhat_tens, 3, symmetrize, simplify = FALSE), along = 3)

  # Returns the estimated precision matrix of the ith observation
  precision <- function(i) {
    Theta <- apply(bhat_symm, c(1, 2), \(b) b %*% c(1, covariates[i, ]))
    diag(Theta) <- 1/varhat
    # omega <- -1 * diag(1/varhat) %*% Theta
    return(Theta)
  }

  # Returns the estimated mean vector of the ith observation
  mean <- function(i) {
    prec <- precision(i)
    mu <- solve(prec, diag(1 / varhat) %*% ghat_mx %*% covariates[i, ])
    return(mu)
  }

  return(list(ghat = ghat_mx,
              bhat = bhat_symm,
              bhat_asym = bhat_tens,
              varhat = varhat,
              lambda = lambdas,
              cv_lambda_idx = cv_lambda_idx,
              precision = precision,
              mean = mean))

  # return(list(ghat = ghat_mx,
  #             bhat = bhat_symm,
  #             bhat_asym = bhat_tens,
  #             varhat = varhat,
  #             precision = precision,
  #             mean = mean))
}