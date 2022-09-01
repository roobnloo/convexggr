library(purrr)

cggr_node <- function(y, responses, covariates,
                      lambda_g, lambda_b_seq, alpha,
                      max_iter = 1000,
                      tol = 1e-5) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(y) == nrow(covariates),
            lambda_g > 0,
            alpha >= 0 && alpha <= 1,
            length(lambda_b_seq) > 0,
            all(lambda_b_seq > 0))

  d <- ncol(responses) + 1
  p <- ncol(covariates)

  # TODO: Once warm starts works, we will sort the weights decreasing.
  # lambda_b_seq <- sort(lambda_b_seq, decreasing = T)

  centered <- center_vars(y, responses, covariates)
  y <- centered$y
  responses <- centered$responses
  covariates <- centered$covariates

  result <- vector(mode = "list", length = length(lambda_b_seq))
  gamma_init <- initialize_gamma(p)
  beta_init <- initialize_beta((p + 1) * (d - 1))
  for (i in seq_along(result)) {
    result[[i]] <- cggr_node_init(y, responses, covariates,
                                  lambda_g, lambda_b_seq[i], alpha,
                                  gamma_init, beta_init, max_iter, tol)
    # TODO: fix the warm starts
    # gamma_init <- result[[i]]$gamma_j
    # beta_init <- result[[i]]$beta_j
  }

  return(result)
}

#' @param y n vector of the response to solve
#' @param reponses n x (d-1) matrix of the other responses
#' @param covariates n x p matrix of covariates
#' @param lambda_g penalty scalar for the mean gamma
#' @param lambda_b penalty scalar for components of beta
#' @param alpha value between 0 and 1 that determines sparse group lasso penalty
cggr_node_init <- function(y, responses, covariates,
                           lambda_g, lambda_b, alpha,
                           gamma_j, beta_j,
                           max_iter, tol) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(y) == nrow(covariates),
            lambda_g > 0,
            lambda_b > 0)

  d <- ncol(responses) + 1
  p <- ncol(covariates)

  r_j <- compute_residual(y, responses, covariates, gamma_j, beta_j)

  obj_val <- numeric(max_iter + 1)
  obj_val[1] <- compute_obj_value(r_j, gamma_j, beta_j, p, d,
                                 lambda_g, lambda_b, alpha)

  for (i in seq_len(max_iter)) {
    result_g <- apply_L1_update(gamma_j, r_j, covariates, lambda_g)
    gamma_j <- result_g$v
    r_j <- result_g$full_resid

    result_b0 <- apply_L1_update(beta_j[seq_len(d-1)], r_j, responses,
                                 lambda_b * alpha)
    beta_j[seq_len(d - 1)] <- result_b0$v
    r_j <- result_b0$full_resid

    for (h in seq_len(p)) {
      bh_idx <- h * (d - 1) + seq_len(d - 1)
      result_bh <- apply_sparsegl_update(beta_j[bh_idx],
                                         r_j, covariates[, h], responses,
                                         lambda_b, alpha)
      beta_j[bh_idx] <- result_bh$bh_j

      ## The code below applies only the LASSO to the groups in beta.
      ## Useful for testing.
      # elt_wise <- apply(responses, 2, function(t) t * covariates[, h])
      # result_bh <- apply_L1_update(beta_j[bh_idx], r_j,
      #                              elt_wise,
      #                              lambda_b * alpha)
      # beta_j[bh_idx] <- result_bh$v
      r_j <- result_bh$full_resid
    }
    obj_val[i + 1] <- compute_obj_value(r_j, gamma_j, beta_j, p, d,
                                       lambda_g, lambda_b, alpha)
    if (is.na(obj_val[i + 1]) | is.nan(obj_val[i + 1])) {
      stop("NaN value encountered when computing L2 loss.
           Perhaps the loss exploded.")
    }
    if (abs(obj_val[i + 1] - obj_val[i]) < tol) {
      obj_val <- obj_val[seq_len(i + 1)]
      break
    }
  }
  if (length(obj_val) == max_iter + 1) {
    warning("Maximum iterations exceeded!")
  }

  sigma_sq <- est_var(r_j, gamma_j, beta_j)

  return(list(gamma_j = gamma_j,
              beta_j = beta_j,
              sigma_sq = sigma_sq,
              obj_val = obj_val,
              resid = r_j))
}

compute_obj_value <- function(r_j, gamma_j, beta_j, p, d,
                              lambda_g, lambda_b, alpha) {
  sum_sq <- sum(r_j^2) / (2 * length(r_j))

  group_lasso_term <-
      map(seq_len(p),
          ~ norm(beta_j[.x * (d - 1) + seq_len(d - 1)], type = "2")) |>
      reduce(`+`)

  sum_sq +
    lambda_g * sum(abs(gamma_j)) +
    alpha * lambda_b * sum(abs(beta_j)) +
    (1 - alpha) * lambda_b * group_lasso_term
}

#' @return scalar estimate of the variance of the error for the jth response
est_var <- function(r_j, gamma_j, beta_j) {
  # threshold the zero coefficients
  tol <- 1e-10
  num_non_zero <- sum(abs(c(gamma_j, beta_j)) < tol)
  sum(r_j^2) / (length(r_j) - num_non_zero)
}

#' Apply the sparse group lasso penalty to encourage both element-wise
#' and group-wise sparsity. The groups are encoded by the interaction between
#' the covariates and the responses.
#' @param bh_j the coefficient vector to be made group-wise and elt-wise sparse.
#' TODO: Most of this function ought to be implemented in C++.
apply_sparsegl_update <- function(bh_j, full_resid, covariate_h, responses,
                                  lambda, alpha) {
  stopifnot(length(full_resid) == nrow(responses),
            length(full_resid) == length(covariate_h))

  if (all(bh_j == 0)) {
    return(list(bh_j = bh_j, full_resid = full_resid))
  }

  n <- nrow(responses)

  grp_intx_mx <- apply(responses, 2, function(x) x * covariate_h)
  grp_partial_resid <- full_resid + grp_intx_mx %*% bh_j
  grp_thresh <- soft_threshold(t(grp_intx_mx) %*% grp_partial_resid / n,
                               alpha * lambda)

  if (sqrt(sum(grp_thresh^2)) <= (1 - alpha) * lambda) {
    return(list(bh_j = 0, full_resid = grp_partial_resid))
  }

  bh_j_norm <- sqrt(sum(bh_j^2))
  bh_j_new <- numeric(length(bh_j))

  for (k in seq_along(bh_j_new)) {
    partial_resid <- full_resid + bh_j[k] * grp_intx_mx[, k]
    inner_prod <- sum(grp_intx_mx[, k] * partial_resid)

    if (abs(inner_prod) <= n * alpha * lambda) {
      bh_j_new[k] <- 0
    } else {
      numerator <- soft_threshold(inner_prod / n, alpha * lambda)
      # it is suggested to multiply the 2nd term in the denom by sqrt(length(bh_j))
      denom <- sum(grp_intx_mx[, k]^2) / n + (1 - alpha) * lambda / bh_j_norm
      bh_j_new[k] <- numerator / denom
    }
  }

  full_resid_new <- full_resid + grp_intx_mx %*% bh_j - grp_intx_mx %*% bh_j_new
  return(list(bh_j = bh_j_new, full_resid = full_resid_new))
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

#' @return n-vector of residuals, where n = length(y)
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
