library(CVXR)
source("cggr_node.R")

cggr_node_cvxr <- function(y, responses, covariates, gamma_j, beta_j,
                                lambda_g, lambda_b, alpha) {
  p <- ncol(covariates)
  d <- ncol(responses) + 1

  # centered <- center_vars(responses, covariates)
  # responses <- centered$responses
  # covariates <- centered$covariates

  loss <- compute_loss(y, responses, covariates, gamma_j, beta_j)

  obj <- loss + cggr_penalty(gamma_j, beta_j, lambda_g, lambda_b, alpha, d, p)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  return(result)
}

compute_loss <- function(y, responses, covariates, gamma_j, beta_j) {
  r_j <- compute_residual(y, responses, covariates, gamma_j, beta_j)
  return(sum_squares(r_j) / (2 * nrow(responses)))
}

cggr_penalty <- function(gamma_j, beta_j, lambda_g, lambda_b, alpha, d, p){
  group_lasso_term <-
      map(seq_len(p), ~ p_norm(beta_j[.x * (d - 1) + seq_len(d - 1)]), 2) |>
        reduce(`+`)

  lambda_g * p_norm(gamma_j, p = 1) +
    alpha * lambda_b * p_norm(beta_j, p = 1) +
    (1 - alpha) * lambda_b * group_lasso_term
}
