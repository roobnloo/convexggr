library(sparsegl)
library(abind)
source("cv_cggr_node.R")
source("node_strategy_sparsegl.R")

#' @param X n x d matrix of responses
#' @param U n x p matrix of covariates
#' @param regmean small ridge penalty on mean vector
#' @param asparse SGL mixture penalty
cggr <- function(X, U, asparse, nodereg_strat = NULL,
                 regmean = 0.01, nlambda = 100) {
  stopifnot(is.matrix(X), is.matrix(U),
            nrow(X) == nrow(U),
            regmean > 0)

  if (is.null(nodereg_strat)) {
    nodereg_strat = node_strategy_sparsegl
  }

  d <- ncol(X)
  p <- ncol(U)
  n <- nrow(X)

  # Tune lambda for all nodes
  lambda <- matrix(nrow = d, ncol = nlambda)
  lambda_min <- vector(length = d)
  for (node in seq_len(d)) {
    result <- cv_cggr_node(node, X, U, asparse, regmean, nodereg_strat,
                           nlambda = nlambda)
    lambda[node,] <- result$lambda
    lambda_min[node] <- result$lambda_min
  }

  # return(list(cv_mse = result$cv_mse,
  #             lambda = lambda,
  #             lambda_min = lambda_min))

  # Initialize mean matrix
  ghat_mx <- matrix(nrow = d, ncol = p)

  # Initialize covariate array
  # Includes the population matrix, hence +1
  bhat_tens <-  array(0, dim = c(d, d, p + 1))

  # Estimated variances
  varhat <- vector(length = d)

  for (node in seq_len(d)) {
    result <- nodereg_strat(node, X, U, lambda_min[node], asparse, regmean)
    rss <- sum(result$resid^2)
    num_nz <- sum(abs(result$beta) > 1e-10) + sum(abs(result$gamma) > 1e-10)
    varhat[node] <- rss  / (n - num_nz)

    bhat_tens[node, -node,] <- -result$beta / varhat[node] # Do we need a scaling factor here? I think so.
    ghat_mx[node,] <- result$gamma
  }

  bhat_symm <- abind(apply(bhat_tens, 3, symmetrize, simplify = F), along = 3)

  # Returns the estimated precision matrix of the ith observation
  precision <- function(i) {
    Theta <- apply(bhat_symm, c(1, 2), \(b) b %*% c(1, U[i, ]))
    diag(Theta) <- 1/varhat
    # omega <- -1 * diag(1/varhat) %*% Theta
    return(Theta)
  }

  # Returns the estimated mean vector of the ith observation
  mean <- function(i) {
    prec <- precision(i)
    mu <- solve(prec) %*% diag(1/varhat) %*% ghat_mx %*% U[i, ]
    return(mu)
  }

  return(list(ghat = ghat_mx,
              bhat = bhat_symm,
              bhat_asym = bhat_tens,
              varhat = varhat,
              precision = precision,
              mean = mean))
}
