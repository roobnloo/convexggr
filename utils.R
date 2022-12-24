#' @return n x p(d-1) matrix representing interactions btw responses and covs
intxmx <- function(responses, covariates) {
  d <- ncol(responses) + 1
  p <- ncol(covariates)
  idx_mat <- as.matrix(expand.grid(seq_len(d - 1), seq_len(p)))

  foo <- function(i) {
    responses[, idx_mat[i, 1]] * covariates[, idx_mat[i, 2]]
  }

  result <- Reduce(
    cbind,
    Map(foo, seq_len(nrow(idx_mat)))
  )
  return(result)
}

#' @return symmetrized version of matrix mx
#' result_ij = result_ji is nonzero iff both mx_ij and mx_ji are nonzero,
#' in which case we choose the smaller value in magnitude.
symmetrize <- function(mx, rule = "and") {
  if (rule == "and") {
    result <- mx * (abs(mx) < t(abs(mx))) + t(mx) * (t(abs(mx)) < abs(mx))
  } else {
    result <- mx * (abs(mx) >= t(abs(mx))) + t(mx) * (t(abs(mx)) >= abs(mx))
  }
  return(result)
}

compute_residual <- function(y, responses, covariates, gamma_j, beta_j) {
  # p <- ncol(responses) + 1
  # wj <- intxmx(responses, covariates)
  # x_gamma_j <- covariates %*% gamma_j # X gamma_j
  # y_b_j0 <- responses %*% beta_j[seq_len(p - 1)] # Y_-j b_j^0
  # wjbetaj0 <- wj %*% beta_j[-seq_len(p - 1)] # W_-j beta_j,-0
  r_beta <- cbind(responses, intxmx(responses, covariates))

  return(y - covariates %*% gamma_j - r_beta %*% beta_j)
}
