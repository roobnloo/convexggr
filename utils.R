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
center_vars <- function(y, responses, covariates) {
  stopifnot(length(y) == nrow(responses),
            nrow(responses) == nrow(covariates))
  n <- length(y)

  scale_n <- function(t) {
    sd(t) * sqrt((n - 1) / n)
  }

  y <- scale(y, scale = FALSE) # center, but do not scale y
  responses <- scale(responses, scale = apply(responses, 2, scale_n))
  covariates <- scale(covariates, scale = apply(covariates, 2, scale_n))

  return(list(y = y,
              responses = responses,
              covariates = covariates))
}
