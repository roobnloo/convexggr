library(CVXR)
source("convexggr.R")

cggr_component_cvxr <- function(y, responses, covariates, gamma_j, beta_j,
                                lambda, alpha) {
  p <- ncol(covariates)
  d <- ncol(responses) + 1
  loss <- compute_loss(y, responses, covariates, gamma_j, beta_j)

  obj <- loss + cggr_penalty(gamma_j, beta_j, lambda, alpha, d, p)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  return(result)
}

compute_loss <- function(y, responses, covariates, gamma_j, beta_j) {
  r_j <- compute_residual(y, responses, covariates, gamma_j, beta_j)
  return(sum_squares(r_j) / (2 * length(r_j)))
}

cggr_penalty <- function(gamma_j, beta_j, lambda, alpha, d, p){
  group_lasso_term <-
      map(seq_len(p), ~ cvxr_norm(beta_j[.x * (d - 1) + seq_len(d - 1)]), 2) |>
        reduce(`+`)

  lambda[1] * cvxr_norm(gamma_j, p = 1) +
    alpha * lambda[2] * cvxr_norm(beta_j, p = 1) +
    (1 - alpha) * lambda[2] * group_lasso_term
}
