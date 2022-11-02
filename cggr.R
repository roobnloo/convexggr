library(sparsegl)
library(abind)
source("cv_cggr_node.R")
source("utils.R")

#' @param X n x d matrix of responses
#' @param U n x p matrix of covariates
#' @param lambda_b_seq vector of d penalties for each response
#' @param alpha_seq values between 0 and 1 that determines sparse group lasso penalty
cggr <- function(X, U, asparse, regmean = 0.01, nlambda = 100) {
  stopifnot(is.matrix(X), is.matrix(U),
            nrow(X) == nrow(U),
            regmean > 0)

  d <- ncol(X)
  p <- ncol(U)
  n <- nrow(X)

  # Tune lambda for all nodes
  lambda <- matrix(nrow = d, ncol = nlambda)
  lambda_min <- vector(length = d)
  for (node in seq_len(d)) {
    result <- cv_cggr_node(node, X, U, asparse, regmean, nlambda = nlambda)
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
    result <- cggr_node(node, X, U, lambda_min[node], asparse, regmean)
    rss <- sum(result$resid^2)
    num_nz <- sum(abs(result$beta) > 1e-10) + sum(abs(result$gamma) > 1e-10)
    varhat[node] <- rss  / (n - num_nz)

    bhat_tens[node, -node,] <- result$beta
#
#     # Index the columns of each beta matrix to fill.
#     # Diagonals are set to zero.
#     beta_mx_idx <- seq_len(d)[-node]
#     for (h in seq_len(p + 1)) {
#       diag(beta_mxs[[h]])[node] <- 0
#       bh_idx <- (h - 1) * (d - 1) + seq_len(d - 1)
#
#       if (reparametrize) {
#         beta_mxs[[h]][node, beta_mx_idx] <- result$beta[bh_idx]
#       } else {
#         beta_mxs[[h]][node, beta_mx_idx] <- -result$beta[bh_idx] / varhat[node]
#       }
#     }

    ghat_mx[node,] <- result$gamma
  }

  bhat_symm <- abind(apply(bhat_tens, 3, symmetrize, simplify = F), along = 3)
  # symmed_beta_mxs <-  map(seq_len(p + 1), ~ matrix(nrow = d, ncol = d))
  # for (h in seq_len(p + 1)) {
  #   symmed_beta_mxs[[h]] <- symmetrize(beta_mxs[[h]])
  # }

  return(list(ghat = ghat_mx,
              bhat = bhat_symm,
              bhat_asym = bhat_tens,
              sigma_sq = varhat))
}
