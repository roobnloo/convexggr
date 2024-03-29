library(modelr)
source("cggr_node.R")

#' @param y n vector of the response to solve
#' @param reponses n x (d-1) matrix of the other responses
#' @param covariates n x p matrix of covariates
#' @param lambda_len length of penalty terms to search
#' @param lambda_min_prop fraction of computed lambda max
#' @param lambda_seq sequence of penalty terms
#' The minimum lambda that yields all zero solutions is calculated by default.
#' Computation can be overwritten by specifying lambda_seq explicitly.
cv_cggr <- function(y, responses, covariates, lambda_g,
                    nfold = 5,
                    alpha_seq = seq(0.1, 1, by = 0.1),
                    lambda_len = 10,
                    lambda_min_prop = 0.1,
                    lambda_seq = NULL) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(y) == nrow(covariates),
            all(alpha_seq > 0) && all(alpha_seq <= 1),
            is.null(lambda_seq) || length(lambda_seq) > 0)

  centered <- center_vars(responses, covariates)
  responses <- centered$responses
  covariates <- centered$covariates


  if(is.null(lambda_seq)) {
    lambda_max <- max(
              sapply(alpha_seq,
                     function(a) min_lambda_zero(y, responses, covariates, a)))
    lambda_seq <- rev(10^seq(log10(lambda_min_prop * lambda_max),
                             log10(lambda_max), length = lambda_len))
  }

  data_df <- as.data.frame(cbind(y, responses, covariates))
  kfold <- crossv_kfold(data_df, nfold)
  fold_mses <- array(NA, dim = c(length(alpha_seq), lambda_len, nfold))
  for (i in seq_along(alpha_seq)) {
    for (j in seq_len(nfold)) {
      fold_mses[i, , j] <- cv_cggr_fold(y, responses, covariates, lambda_g,
                                        lambda_seq, alpha_seq[i],
                                        kfold$train[[j]]$idx,
                                        kfold$test[[j]]$idx)
    }
  }

  cv_error <- apply(fold_mses, c(1, 2), mean)
  opt_idx <- which(cv_error == min(cv_error), arr.ind = T)[1, ]
  return(list(cv_error = cv_error,
              alpha_seq = alpha_seq,
              lambda_seq = lambda_seq,
              opt = list(alpha = alpha_seq[opt_idx[1]],
                         lambda = lambda_seq[opt_idx[2]])))
}

#' @return Numeric vector of length(lambda_seq) giving test mse for a provided
#' train/test split
cv_cggr_fold <- function(y, responses, covariates,
                         lambda_g, lambda_seq, alpha,
                         train_idx, test_idx) {
  y_train <- y[train_idx]
  responses_train <- responses[train_idx, ]
  covariates_train <- covariates[train_idx, ]

  n_train <- length(train_idx)
  sd_n <- function(t) {
    sd(t) * sqrt((n_train - 1) / n_train)
  }

  # List of length(lambda_b_seq)
  fit_seq <- cggr_node(y_train, responses_train, covariates_train,
                       lambda_g, lambda_seq, alpha)

  predict_mse <- numeric(length(lambda_seq))

  # For each result in the list, predict on the test sample
  for (i in seq_along(fit_seq)) {
    # Center and scale the test covariates
    responses_test <- responses[test_idx, ] |>
                      sweep(2, apply(responses_train, 2, mean), "-") |>
                      sweep(2, apply(responses_train, 2, sd_n), "/")
    covariates_test <- covariates[test_idx, ] |>
                       sweep(2, apply(covariates_train, 2, mean), "-") |>
                       sweep(2, apply(covariates_train, 2, sd_n), "/")

    resid <- compute_residual(y[test_idx], responses_test, covariates_test,
                              fit_seq[[i]]$gamma_j, fit_seq[[i]]$beta_j)

    predict_mse[i] <- mean(resid^2)
  }

  return(predict_mse)
}

min_lambda_zero <- function(y, responses, covariates, alpha) {
  d <- ncol(responses) + 1
  p <- ncol(covariates)

  intx <- interaction_mx(responses, covariates)

  lambda_grp <- numeric(p)
  for (h in seq_len(p)) {
    bh_idx <- (h - 1) * (d - 1) + seq_len(d - 1)
    quad_grad <- t(intx[, bh_idx]) %*% y / nrow(responses)
    upper <- norm(quad_grad, "2") / alpha

    piece_quad <- function(t) {
      sum(soft_threshold(quad_grad, t * alpha)^2) - (1 - alpha)^2 * t^2
    }

    lambda_grp[h] <- uniroot(piece_quad, c(0, upper))$root
  }

  return(max(lambda_grp))
}
