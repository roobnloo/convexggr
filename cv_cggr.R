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
cv_cggr <- function(y, responses, covariates, lambda_g, alpha,
                    nfold = 5,
                    lambda_len = 20,
                    lambda_min_prop = 0.1,
                    lambda_seq = NULL) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(y) == nrow(covariates),
            is.null(lambda_seq) || length(lambda_seq) > 0)

  centered <- center_vars(y, responses, covariates)
  y <- centered$y
  responses <- centered$responses
  covariates <- centered$covariates

  if(is.null(lambda_seq)) {
    lambda_max <- min_lambda_zero(y, responses, covariates, alpha)
    lambda_seq <- rev(10^seq(log10(lambda_min_prop * lambda_max),
                             log10(lambda_max), length = lambda_len))
  }

  data_df <- as.data.frame(cbind(y, responses, covariates))
  kfold <- crossv_kfold(data_df, nfold)
  fold_mses <- matrix(NA, nrow = nfold, ncol = length(lambda_seq))
  for (i in seq_len(nfold)) {
    fold_mses[i, ] <- cv_cggr_fold(y, responses, covariates,
                                   lambda_g, lambda_seq, alpha,
                                   kfold$train[[i]]$idx, kfold$test[[i]]$idx)
  }

  result <- apply(fold_mses, 2, mean)
  return(list(result = result,
              lambda_seq = lambda_seq,
              lambda_opt = lambda_seq[which.min(result)]))
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
    # Center the test response and center and scale the test covariates
    y_test <- y[test_idx] - mean(y_train)
    responses_test <- responses[test_idx, ] |>
                      sweep(2, apply(responses_train, 2, mean), "-") |>
                      sweep(2, apply(responses_train, 2, sd_n), "/")
    covariates_test <- covariates[test_idx, ] |>
                       sweep(2, apply(covariates_train, 2, mean), "-") |>
                       sweep(2, apply(covariates_train, 2, sd_n), "/")

    resid <- compute_residual(y_test, responses_test, covariates_test,
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
  for (h in seq_along(p)) {
    bh_idx <- h * (d - 1) + seq_len(d - 1)
    quad_grad <- t(intx[, bh_idx]) %*% y / nrow(responses)
    upper <- norm(quad_grad, "2") / alpha

    piece_quad <- function(t) {
      sum(soft_threshold(quad_grad, t * alpha)^2) - (1 - alpha)^2 * t^2
    }

    lambda_grp[h] <- uniroot(piece_quad, c(0, upper))$root
  }

  return(max(lambda_grp))
}
