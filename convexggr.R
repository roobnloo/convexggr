library(tibble)
library(purrr)

#' @param reponses n x d matrix of responses
#' @param covariates n x p matrix of covariates
#' @param lambda 2 vector of penalty terms
#' @param alpha value between 0 and 1 that determines sparse group lasso penalty
convex_ggr <- function(responses, covariates, lambda, alpha = 0.5,
                       gamma_init = NULL,
                       beta_init = NULL,
                       max_iter = 1000,
                       tol = 1e-5) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(lambda) == 2)

  result <- vector("list", length = ncol(responses))

  # Initialize mean matrix, gamma
  gamma_mx <- matrix(nrow = ncol(responses), ncol = ncol(covariates))

  # Initialize covariate coef mxs, beta.
  # Includes the population matrix, hence +1
  beta_mxs <-  map(seq_len(ncol(covariates) + 1),
                   ~ matrix(nrow = ncol(responses), ncol = ncol(responses)))

  for (i in seq_along(result)) {
    result[[i]] <- convex_ggr_component(responses[, i], responses[, -i],
                                        covariates, lambda, alpha,
                                        gamma_init, beta_init, max_iter, tol)
    gamma_mx[i, ] <- result[[i]]$gamma_j
    beta <- result[[i]]$beta_j

    #(beta_mxs[[1]])[1, ] <- beta[1:]
  }

  return(result)
}

#' @param y n vector of the response to solve
#' @param reponses n x (d-1) matrix of the other responses
#' @param covariates n x p matrix of covariates
#' @param lambda 2 vector of penalty terms
#' @param alpha value between 0 and 1 that determines sparse group lasso penalty
convex_ggr_component <- function(y, responses, covariates, lambda, alpha = 0.5,
                       gamma_init = NULL,
                       beta_init = NULL,
                       max_iter = 1000,
                       tol = 1e-5) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(y) == nrow(covariates),
            length(lambda) == 2)

  n <- nrow(responses)
  d <- ncol(responses) + 1
  p <- ncol(covariates)

  centered <- center_vars(y, responses, covariates)
  y <- centered$y
  responses <- centered$responses
  covariates <- centered$covariates

  gamma_j <- gamma_init
  beta_j <- beta_init

  if (is.null(gamma_j)) {
    gamma_j <- initialize_gamma(p)
  }
  if (is.null(beta_j)) {
    beta_j <- initialize_beta((p + 1) * (d - 1))
  }

  r_j <- compute_residual(y, responses, covariates, gamma_j, beta_j)

  losses <- numeric(max_iter + 1)
  losses[1] <- sq_error_loss(r_j)

  for (i in seq_len(max_iter)) {
    result_g <- apply_L1_update(gamma_j, r_j, covariates, lambda[1])
    gamma_j <- result_g$v
    r_j <- result_g$full_resid

    result_b0 <- apply_L1_update(beta_j[seq_len(d-1)], r_j, responses,
                                 lambda[2] * alpha)
    beta_j[seq_len(d - 1)] <- result_b0$v
    r_j <- result_b0$full_resid

    for (h in seq_len(p)) {
      bh_idx <- h * (d - 1) + seq_len(d - 1)
      result_bh <- apply_sparsegl_update(beta_j[bh_idx],
                                         r_j, covariates[, h], responses,
                                         lambda[2], alpha)
      beta_j[bh_idx] <- result_bh$bh_j
      r_j <- result_bh$full_resid
    }

    losses[i + 1] <- sq_error_loss(r_j)
    if (is.na(losses[i + 1]) | is.nan(losses[i + 1])) {
      stop("NaN value encountered when computing L2 loss.
           Perhaps the loss exploded.")
    }
    if (abs(losses[i + 1] - losses[i]) < tol) {
      losses <- losses[seq_len(i + 1)]
      break
    }
  }
  if (length(losses) == max_iter + 1) {
    warning("Maximum iterations exceeded!")
  }
  return(list(gamma_j = gamma_j,
              beta_j = beta_j,
              losses = losses,
              resid = r_j))
}

#' Apply the sparse group lasso penalty to encourage both element-wise
#' and group-wise sparsity. The groups are encoded by the interaction between
#' the covariates and the responses.
#' @param bh_j the coefficient vector to be made group-wise and elt-wise sparse.
apply_sparsegl_update <- function(bh_j, full_resid, covariate_h, responses,
                                  lambda0, alpha) {
  stopifnot(length(full_resid) == nrow(responses),
            length(full_resid) == length(covariate_h))
  n <- nrow(responses)
  bh_j_new <- numeric(length(bh_j))
  full_resid_new <- full_resid
  bh_j_norm <- norm(bh_j, "2")

  for (k in seq_len(ncol(responses))) {
    element_wise <- responses[, k] * covariate_h
    partial_resid <- full_resid + bh_j[k] * element_wise

    numerator <- sum(element_wise * partial_resid) / n |>
                   soft_threshold(alpha * lambda0)

    # it is suggested to multiply the second term in the denom by sqrt(length(bh_j))
    denom <- sum(element_wise^2) / n +
              (1 - alpha) * lambda0 / bh_j_norm

    bh_j_new[k] <- numerator / denom

    full_resid <- partial_resid - bh_j[k] * element_wise
    full_resid_new <- full_resid_new + bh_j[k] * element_wise -
                      bh_j_new[k] * element_wise
  }

  list(bh_j = bh_j_new, full_resid = full_resid_new)
}

#' Applies the coordinate-wise update for the L1 penalty
#' This encourages element-wise sparsity on the vector.
#' @param v the coefficient vector to be made sparse.
#' @param lambda the penalty coefficient. Larger values mean more sparsity.
#' @param design_mx design matrix which multiplies with coefficient vector v
apply_L1_update <- function(v, full_resid, design_mx, lambda) {
  stopifnot(length(full_resid) == nrow(design_mx),
            length(v) == ncol(design_mx))

  for (k in seq_along(v)) {
    partial_resid <- full_resid + v[k] * design_mx[, k]
    v[k] <- soft_threshold(v[k] + sum(design_mx[, k] * full_resid) / nrow(design_mx),
                           lambda)
    full_resid <- partial_resid - v[k] * design_mx[, k]
  }

  list(v = v, full_resid = full_resid)
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
  stopifnot(length(lambda) == 1)
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

reshape <- function(beta_j) {

}
