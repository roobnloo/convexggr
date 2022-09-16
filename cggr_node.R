library(purrr)
source("utils.R")

cggr_node <- function(y, responses, covariates,
                      lambda_g, lambda_b_seq, alpha,
                      max_iter = 1000,
                      tol = 1e-5) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(y) == nrow(covariates),
            lambda_g > 0,
            length(alpha) == 1,
            alpha >= 0 && alpha <= 1,
            length(lambda_b_seq) > 0,
            all(lambda_b_seq > 0))

  d <- ncol(responses) + 1
  p <- ncol(covariates)

  centered <- center_vars(responses, covariates)
  responses <- centered$responses
  covariates <- centered$covariates

  result <- vector(mode = "list", length = length(lambda_b_seq))
  gamma_init <- rep(0, p)
  beta_init <- rep(0, (p + 1) * (d - 1))
  for (i in seq_along(result)) {
    result[[i]] <- cggr_node_init(y, responses, covariates,
                                  lambda_g, lambda_b_seq[i], alpha,
                                  gamma_init, beta_init, max_iter, tol)

    # Initialize with previous results (warm starts)
    gamma_init <- result[[i]]$gamma_j
    beta_init <- result[[i]]$beta_j
  }

  return(result)
}

#' @param y n vector of the response to solve
#' @param reponses n x (d-1) matrix of the other responses
#' @param covariates n x p matrix of covariates
#' @param lambda_g penalty scalar for the mean gamma
#' @param lambda_b penalty scalar for components of beta
#' @param alpha value between 0 and 1 that determines sparse group lasso penalty
#' @param gamma_j initial gamma_j
#' @param beta initial beta_j
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
  group_lasso_term <-
      map(seq_len(p),
          ~ norm(beta_j[.x * (d - 1) + seq_len(d - 1)], type = "2")) |>
      reduce(`+`)

  quad_loss(r_j) +
    lambda_g * sum(abs(gamma_j)) +
    alpha * lambda_b * sum(abs(beta_j)) +
    (1 - alpha) * lambda_b * group_lasso_term
}

quad_loss <- function(r_j) {
  sum(r_j^2) / (2 * length(r_j))
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
#' @returns updated coefficient vector and residual vector
#' TODO: Most of this function ought to be implemented in C++.
apply_sparsegl_update <- function(bh_j, full_resid, covariate_h, responses,
                                  lambda, alpha) {
  stopifnot(length(full_resid) == nrow(responses),
            length(full_resid) == length(covariate_h))

  n <- nrow(responses)
  grp_intx_mx <- apply(responses, 2, function(x) x * covariate_h)
  grp_partial_resid <- full_resid + grp_intx_mx %*% bh_j
  grp_thresh <- soft_threshold(t(grp_intx_mx) %*% grp_partial_resid / n,
                               alpha * lambda)

  # If this subgradient condition holds, the entire group should be zero
  if (sqrt(sum(grp_thresh^2)) <= (1 - alpha) * lambda) {
    return(list(bh_j = 0, full_resid = grp_partial_resid))
  }

  step_size <- 1
  bh_j_new <- bh_j
  theta_new <- bh_j
  maxit <- 100
  obj <- Inf
  for (l in seq_len(maxit)) {
    # TODO: clean this up
    obj_new <- compute_obj_value(grp_partial_resid - grp_intx_mx %*% bh_j_new,
                                  0, c(rep(0, length(bh_j_new)), bh_j_new),
                                  1, length(bh_j_new) + 1, 0, lambda, alpha)
    if (abs(obj_new - obj) < 1e-5) {
      break
    }
    if (l == maxit) {
      warning(paste("Inner loop within SGL descent exceeded max", maxit, "iterations"))
    }

    obj <- obj_new
    theta_old <- theta_new
    grp_fit <- grp_intx_mx %*% bh_j_new
    grad <- -1 * t(grp_intx_mx) %*% (grp_partial_resid - grp_fit) / n

    # Optimize the step size
    quad_loss_old <- quad_loss(grp_partial_resid - grp_fit)
    repeat {
      theta_new <- sgl_gd_update(bh_j_new, step_size, grad, alpha, lambda)
      delta <- theta_new - bh_j_new
      rhs <- quad_loss_old + sum(grad * delta) + sum(delta^2) / (2 * step_size)
      lhs <- quad_loss(grp_partial_resid - grp_intx_mx %*% theta_new)
      if (lhs <= rhs) {
        break
      }
      step_size <- step_size * 0.8
    }

    # Nesterov momentum step
    bh_j_new <- theta_old + (l / (l + 3)) * (theta_new - theta_old)
  }

  full_resid_new <- grp_partial_resid - grp_intx_mx %*% bh_j_new
  return(list(bh_j = as.numeric(bh_j_new), full_resid = full_resid_new))
}

sgl_gd_update <- function(center, step_size, grad, alpha, lambda) {
    thresh_step <- soft_threshold(center - step_size * grad,
                                  step_size * alpha * lambda)
    ss <- sum(thresh_step^2)
    denom <- ifelse(ss == 0, 1, sqrt(ss))
    max_term <- max(0, 1 - step_size * (1 - alpha) * lambda / denom)
    return(max_term * thresh_step)
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

# Vectorized soft-threshold function
soft_threshold <- function(x, lambda) {
  stopifnot(length(lambda) == 1)
  sign(x) * pmax(abs(x) - lambda, 0)
}
