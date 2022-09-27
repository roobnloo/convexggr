library(purrr)

#' @return n x p(d-1) matrix representing interactions btw responses and covs
interaction_mx <- function(responses, covariates) {
  d <- ncol(responses) + 1
  p <- ncol(covariates)
  idx_mat <- as.matrix(expand.grid(seq_len(d - 1), seq_len(p)))

  foo <- function(i) {
    responses[, idx_mat[i, 1]] * covariates[, idx_mat[i, 2]]
  }

  map(seq_len(nrow(idx_mat)), foo) |>
    reduce(cbind)
}

#' @return variables with mean zero and sum-of-squares equal to nrow
center_vars <- function(responses, covariates) {
  stopifnot(nrow(responses) == nrow(covariates))
  n <- nrow(responses)

  scale_n <- function(t) {
    sd(t) * sqrt((n - 1) / n)
  }

  responses <- scale(responses, scale = apply(responses, 2, scale_n))
  covariates <- scale(covariates, scale = apply(covariates, 2, scale_n))

  return(list(responses = responses,
              covariates = covariates))
}

#' @return n-vector of residuals, where n = length(y)
compute_residual <- function(y, responses, covariates, gamma_j, beta_j) {
  d <- ncol(responses) + 1
  W_j <- interaction_mx(responses, covariates)
  x_gamma_j <- covariates %*% gamma_j # X gamma_j
  y_b_j0 <- responses %*% beta_j[seq_len(d-1)] # Y_-j b_j^0
  W_jbeta_j0 <- W_j %*% beta_j[-seq_len(d-1)] # W_-j beta_j,-0

  return(y - x_gamma_j - y_b_j0 - W_jbeta_j0)
}

#' @return symmetrized version of matrix mx
#' result_ij = result_ji is nonzero iff both mx_ij and mx_ji are nonzero,
#' in which case we choose the smaller value in magnitude.
symmetrize <- function(mx) {
  ut <- mx[upper.tri(mx)]
  lt <- t(mx)[upper.tri(mx)]

  symmed <- unlist(map2(ut, lt, symm_help))
  mx[upper.tri(mx)] <- symmed
  for(i in seq_len(nrow(mx))) {
    for(j in seq_len(i - 1)){
      mx[i, j] <- mx[j, i]
    }
  }
  return(mx)
}

symm_help <- function(x, y) {
  if (isTRUE(all.equal(x, 0)) & isTRUE(all.equal(y, 0))) {
    return(0)
  }
  ifelse(abs(x) < abs(y), x, y)
}
