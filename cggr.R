library(abind)
library(Rcpp)
source("utils.R")
source("cv_cggr_node.R")
# source("node_strategy_sparsegl.R")
sourceCpp("nodewiseRegression.cpp")

#' @param responses n x d matrix of responses
#' @param U n x p matrix of covariates
#' @param nregmean
#' @param nlambda
#' @param asparse SGL mixture penalty
cggr <- function(responses, covariates, asparse,
                 nregmean = 10, nlambda = 100, lambdafactor = 1e-4,
                 maxit = 3e6, tol = 1e-8, nfolds = 5,
                 verbose = FALSE, parallel = TRUE) {

  stopifnot(is.matrix(responses), is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            nregmean > 0)

  p <- ncol(responses)
  q <- ncol(covariates)
  n <- nrow(responses)
  bveclength <- (p - 1) * (q + 1)

  # Initial run to get lambdas for each response
  lambdas <- matrix(nrow = nlambda, ncol = p)
  regmeanpaths <- matrix(nrow = nregmean, ncol = p)
  beta <- array(dim = c(bveclength, nlambda, nregmean, p))
  gamma <- array(dim = c(q, nlambda, nregmean, p))
  varhat <- array(dim = c(nlambda, nregmean, p))
  resid <- array(dim = c(n, nlambda, nregmean, p))
  objval <- array(dim = c(nlambda, nregmean, p))

  nodewise <- function(node) {
    # if (verbose)
    #   print(paste("Starting initial run for node", node))
    nodereg <- nodewiseRegression(
      responses[, node], responses[, -node], covariates, asparse,
      nregmean = nregmean, nlambda = nlambda, lambdaFactor = lambdafactor,
      maxit = maxit, tol = tol)
    if (verbose)
      print(paste("Finished initial run for node", node))
    return(list(
      lambdas = nodereg["lambdas"][[1]],
      regmeanpath = nodereg["regmeans"][[1]],
      beta = nodereg["beta"][[1]],
      gamma = nodereg["gamma"][[1]],
      varhat = nodereg["varhat"][[1]],
      resid = nodereg["resid"][[1]],
      objval = nodereg["objval"][[1]]
    ))
  }

  if (verbose) {
    cat("Begin initial run...\n")
    tictoc::tic()
  }
  if (parallel) {
    reg_result <- parallel::mclapply(seq_len(p), nodewise, mc.cores = 10L)
  } else {
    reg_result <- lapply(seq_len(p), nodewise)
  }

  for (node in seq_len(p)) {
    lambdas[, node] <- reg_result[[node]]$lambdas
    regmeanpaths[, node] <- reg_result[[node]]$regmeanpath
    beta[, , , node] <- reg_result[[node]]$beta
    gamma[, , , node] <- reg_result[[node]]$gamma
    varhat[, , node] <- reg_result[[node]]$varhat
    resid[, , , node] <- reg_result[[node]]$resid
    objval[, , node] <- reg_result[[node]]$objval
  }

  if (verbose) {
    print("Finished initial run")
    tictoc::toc()
    print("Begin cross-validation...")
    tictoc::tic()
  }

  cv_lambda_idx <- vector(length = p)
  cv_regmean_idx <- vector(length = p)
  cv_mse <- array(dim = c(p, nlambda, nregmean))

  cv_node <- function(node) {
    cv_result <- cv_cggr_node(
      node, responses, covariates, asparse,
      lambdas[, node], regmeanpaths[, node],
      maxit, tol, nfolds)
    if (verbose)
      print(paste("Done cross-validating node", node))
    return(cv_result)
  }

  if (parallel) {
    cv_results <- parallel::mclapply(seq_len(p), cv_node, mc.cores = 10L)
  } else {
    cv_results <- lapply(seq_len(p), cv_node)
  }

  # if (any(sapply(cv_results, inherits, what = "try-error"))) {
  #   browser()
  # }

  for (node in seq_len(p)) {
    cv_mse[node, , ] <- cv_results[[node]]
    minind <- arrayInd(which.min(cv_results[[node]]), dim(cv_results[[node]]))
    cv_lambda_idx[node] <- minind[1]
    cv_regmean_idx[node] <- minind[2]
  }

  if (verbose) {
    print("Finished cross validating all nodes")
    tictoc::toc()
  }

  ghat_select <- cbind(rep(seq(q), times = p),
                       rep(cv_lambda_idx, each = q),
                       rep(cv_regmean_idx, each = q),
                       rep(seq(p), each = q))
  ghat_mx <- t(matrix(gamma[ghat_select], nrow = q, ncol = p))

  varhat <- varhat[cbind(cv_lambda_idx, cv_regmean_idx, seq(p))]

  bhat_select <- cbind(rep(seq(bveclength), times = p),
                       rep(cv_lambda_idx, each = bveclength),
                       rep(cv_regmean_idx, each = bveclength),
                       rep(seq(p), each = bveclength))
  bhat_mx <- matrix(beta[bhat_select], nrow = bveclength, ncol = p)
  bhat_tens <-  array(0, dim = c(p, p, q + 1))
  for (i in seq_len(p)) {
    bhat_tens[i, -i, ] <- bhat_mx[, i] / varhat[i]
  }

  bhat_symm <- abind(
    apply(bhat_tens, 3, symmetrize, simplify = FALSE), along = 3)

  # Returns the estimated precision matrix of the ith observation
  precision <- function(i) {
    Theta <- apply(bhat_symm, c(1, 2), \(b) b %*% c(1, covariates[i, ]))
    omega <- -1 * Theta
    diag(omega) <- 1 / varhat
    return(omega)
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
              regmeanpaths = regmeanpaths,
              cv_lambda_idx = cv_lambda_idx,
              cv_regmean_idx = cv_regmean_idx,
              precision = precision,
              mean = mean))
}