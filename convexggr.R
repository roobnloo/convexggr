library(tibble)
library(purrr)

#' @param y n vector of the response to solve
#' @param reponses n x (d-1) matrix of the other responses
#' @param covariates n x p matrix of covariates
convex_ggr <- function(y, responses, covariates, max_iter = 1000) {
  stopifnot(is_matrix(responses),
            is_matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(y) == nrow(covariates))
  n <- nrow(responses)
  d <- ncol(responses) + 1
  p <- ncol(covariates)
  W_j <- interaction_mx(responses, covariates)

  gamma_j <- initialize_gamma(p)
  beta_j <- initialize_beta((p + 1) * (d - 1))

  losses <- numeric(1000)
  for (i in seq_len(max_iter)) {
    loss[i] <- sq_error_loss(y, responses, covariates, gamma_j, beta_j)

  }
}

#' @return n x p(d-1) matrix representing interactions btw responses and covs
interaction_mx <- function(responses, covariates) {
  d <- ncol(responses) + 1
  p <- ncol(covariates)
  idx_mat <- as.matrix(expand.grid(seq_len(d - 1), seq_len(p)))

  map(seq_len(nrow(idx_mat)), ~ responses[, idx_mat[.x, 1]] *
                                covariates[, idx_mat[.x, 2]]) |>
    reduce(cbind)
}

#' @return squared L2 loss given parameters gamma_j and beta_j
sq_error_loss <- function(y, responses, covariates, gamma_j, beta_j) {
  d <- ncol(responses) + 1
  W_j <- interaction_mx(responses, covariates)
  x_gamma_j <- covariates %*% gamma_j # X gamma_j
  y_b_j0 <- responses %*% beta_j[seq_len(d-1)] # Y_-j b_j^0
  W_jbeta_j0 <- W_j %*% beta_j[-seq_len(d-1)] # W_-j beta_j,-0

  norm(y - x_gamma_j - y_b_j0 - W_jbeta_j0, type = "f")^2 / (2 * nrow(responses))
}

# Vectorized soft-threshold function
soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

#' @return initialized p-vector gamma_vec
initialize_gamma <- function(p) {
  runif(p)
}

#' @return initialized q-vector beta_vec
initialize_beta <- function(q) {
  runif(q)
}
