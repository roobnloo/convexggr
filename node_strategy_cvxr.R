library(CVXR)
source("utils.R")

node_strategy_cvxr <- function(node, X, U, lambda, asparse, regmean,
                               initbeta, initgamma,
                               maxit = 1000, tol = 1e-5) {

  p <- ncol(U)
  d <- ncol(X)
  y <- scale(X[, node], scale = F)
  X <- X[, -node]

  gamma <- initgamma
  beta <- initbeta
  loss <- compute_loss(y, X, U, gamma, beta)

  obj <- loss + cggr_penalty(gamma, beta, regmean, lambda, asparse, d, p)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  return(result)
}

compute_loss <- function(y, responses, covariates, gamma_j, beta_j) {
  r_j <- compute_residual(y, responses, covariates, gamma_j, beta_j)
  return(sum_squares(r_j) / (2 * nrow(responses)))
}

cggr_penalty <- function(gamma_j, beta_j, lambda_g, lambda_b, alpha, d, p) {
  group_lasso_term <-
      map(seq_len(p), ~ p_norm(beta_j[.x * (d - 1) + seq_len(d - 1)]), 2) |>
        reduce(`+`)

  lambda_g * sum_squares(gamma_j) +
    alpha * lambda_b * p_norm(beta_j, p = 1) +
    (1 - alpha) * lambda_b * group_lasso_term
}
