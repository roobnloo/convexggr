library(tibble)
library(purrr)

#' @param y n vector of the response to solve
#' @param reponses n x (d-1) matrix of the other responses
#' @param covariates n x p matrix of covariates
#' @param lambda 2 vector of penalty terms
#' @param alpha value between 0 and 1 that determines sparse group lasso penalty
convex_ggr <- function(y, responses, covariates, lambda, alpha = 0.5,
                       gamma_init = NULL,
                       beta_init = NULL,
                       max_iter = 1000,
                       tol = 1e-5) {
  # stopifnot(is.matrix(responses),
  #           is.matrix(covariates),
  #           nrow(responses) == nrow(covariates),
  #           length(y) == nrow(covariates),
  #           length(lambda == 2))
  n <- nrow(responses)
  d <- ncol(responses) + 1
  p <- ncol(covariates)
  # W_j <- interaction_mx(responses, covariates)

  gamma_j <- gamma_init
  beta_j <- beta_init

  if (is.null(gamma_j)) {
    gamma_j <- initialize_gamma(p)
  }
  if (is.null(beta_j)) {
    beta_j <- initialize_beta((p + 1) * (d - 1))
  }

  r_j <- compute_residual(y, responses, covariates, gamma_j, beta_j)

  losses <- numeric(max_iter)
  losses[1] <- sq_error_loss(r_j)

  for (i in seq_len(max_iter - 1) + 1) {
    result_g <- update_gamma(gamma_j, r_j, covariates, lambda[1])
    gamma_j <- result_g$gamma_j
    r_j <- result_g$full_resid

    result_b0 <- update_b0(beta_j[seq_len(d-1)], r_j, responses,
                           lambda[2], alpha)
    beta_j[seq_len(d-1)] <- result_b0$b0_j
    r_j <- result_b0$full_resid

    for (h in seq_len(p)) {
      result_bh <- update_bh(beta_j[h * (d - 1) + seq_len(d-1)], r_j,
                             covariates[, h], responses, lambda[2], alpha)
      beta_j[h * (d - 1) + seq_len(d-1)] <- result_bh$bh_j
      r_j <- result_bh$full_resid
    }
    losses[i] <- sq_error_loss(r_j)
    if (is.na(losses[i])) {
      browser()
    }
    if (abs(losses[i] - losses[i - 1]) < tol) {
      return(list(gamma_j = gamma_j,
                  beta_j = beta_j,
                  losses = losses,
                  resid = r_j))
    }
  }
}

update_bh <- function(bh_j, full_resid, covariate_h, responses, lambda0, alpha) {
  stopifnot(length(full_resid) == nrow(responses),
            length(full_resid) == length(covariate_h))

  n <- nrow(responses)
  browser()

  for (k in seq_len(ncol(responses))) {
    element_wise <- responses[, k] * covariate_h
    partial_resid <- full_resid + bh_j[k] * element_wise

    numerator <- sum(element_wise * partial_resid) |>
                   soft_threshold(alpha * lambda0)
    denom <- norm(element_wise, "2")^2 / n +
              (1 - alpha) * lambda0 / norm(bh_j, "2")

    bh_j[k] <- numerator / denom

    full_resid <- partial_resid - bh_j[k] * element_wise
  }

  list(bh_j = bh_j, full_resid = full_resid)
}

update_b0 <- function(b0_j, full_resid, responses, lambda0, alpha) {
  stopifnot(length(full_resid) == nrow(responses),
            length(b0_j) == ncol(responses))

  n <- nrow(responses)

  for (k in seq_len(ncol(responses))) {
    partial_resid <- full_resid + b0_j[k] * responses[, k]
    b0_j[k] <- b0_j[k] + sum(responses[, k] * partial_resid) / n |>
                  soft_threshold(alpha * lambda0)
    full_resid <- partial_resid - b0_j[k] * responses[, k]
  }

  list(b0_j = b0_j, full_resid = full_resid)
}

update_gamma <- function(gamma_j, full_resid, covariates, lambda1) {
  stopifnot(length(full_resid) == nrow(covariates),
            length(gamma_j) == ncol(covariates))

  n <- nrow(covariates)

  for (k in seq_len(ncol(covariates))) {
    partial_resid <- full_resid + gamma_j[k] * covariates[, k]
    gamma_j[k] <- soft_threshold(gamma_j[k] +
                                   sum(covariates[, k] * partial_resid) / n,
                                 lambda1)
    full_resid <- partial_resid - gamma_j[k] * covariates[, k]
  }

  list(gamma_j = gamma_j, full_resid = full_resid)
}

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

#' @return squared L2 loss given parameters gamma_j and beta_j
sq_error_loss <- function(full_residual) {
  stopifnot(is.atomic(full_residual))
  norm(full_residual, type = "2")^2 / (2 * length(full_residual))
}

compute_residual <- function(y, responses, covariates, gamma_j, beta_j) {
  d <- ncol(responses) + 1
  W_j <- interaction_mx(responses, covariates)
  x_gamma_j <- covariates %*% gamma_j # X gamma_j
  y_b_j0 <- responses %*% beta_j[seq_len(d-1)] # Y_-j b_j^0
  W_jbeta_j0 <- W_j %*% beta_j[-seq_len(d-1)] # W_-j beta_j,-0

  return(y - x_gamma_j - y_b_j0 - W_jbeta_j0)
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
