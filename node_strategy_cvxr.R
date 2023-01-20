library(CVXR)
source("utils.R")

node_strategy_cvxr <- function(y, responses, covariates,
			       lambda, asparse, regmean,
			       initbeta, initgamma,
			       maxit = 1000, tol = 1e-5) {

  q <- ncol(covariates)
  p <- ncol(responses) + 1

  gamma <- initgamma
  beta <- initbeta
  loss <- compute_loss(y, responses, covariates, gamma, beta)

  obj <- loss + cggr_penalty(gamma, beta, regmean, lambda, asparse, p, q)
  prob <- Problem(Minimize(obj))
  result <- CVXR::solve(
    prob, reltol = tol, abstol = tol, num_iter = maxit,
    verbose = TRUE)
  return(result)
}

compute_loss <- function(y, responses, covariates, gamma_j, beta_j) {
  r_j <- compute_residual(y, responses, covariates, gamma_j, beta_j)
  return(sum_squares(r_j) / (2 * nrow(responses)))
}

cggr_penalty <- function(gamma_j, beta_j, lambda_g, lambda_b, alpha, p, q) {
  group_lasso_term <-
    Map(\(x) p_norm(beta_j[x * (p - 1) + seq_len(p - 1)], p = 2), seq_len(q))
  group_lasso_term <- Reduce(`+`, group_lasso_term)

  # gl <- sum(sapply(seq_len(q),
  #                  \(x) p_norm(beta_j[x * (p-1) + seq_len(p-1)], p = 2)
  #          ))

  lambda_g * sum_squares(gamma_j) +
    alpha * lambda_b * p_norm(beta_j, p = 1) +
    (1 - alpha) * lambda_b * group_lasso_term
}
