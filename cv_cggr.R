library(modelr)
source("cggr_node.R")

#' @param y n vector of the response to solve
#' @param reponses n x (d-1) matrix of the other responses
#' @param covariates n x p matrix of covariates
#' @param lambda_seq sequence of penalty terms
cv_cggr <- function(y, responses, covariates, lambda_seq, nfold = 5) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            length(y) == nrow(covariates),
            length(lambda_seq) > 0)

  data_df <- as.data.frame(cbind(y, responses, covariates))
  kfold <- crossv_kfold(data_df, nfold)

  # gamma_init <- rep(0, p)
  # beta_init <- rep(0, d - 1)

  # cv_mses <- sapply(lambda_seq,
  #                function (t) cv_cggr_pen(y, responses, covariates, kfold, t))

  fold_mses <- matrix(NA, nrow = nfold, ncol = length(lambda_seq))
  for (i in seq_len(nfold)) {
    fold_mses[i, ] <- cv_cggr_fold(y, responses, covariates, lambda_seq,
                                   kfold$train[[i]]$idx, kfold$test[[i]]$idx)
  }

  result <- apply(fold_mses, 2, mean)
  return(result)
}

#' @return Numeric vector of length(lambda_seq) giving test mse for a provided
#' train/test split
cv_cggr_fold <- function(y, responses, covariates, lambda_seq,
                         train_idx, test_idx,
                         lambda_g = 0.1, alpha = 0.5) {
  y_train <- y[train_idx]
  responses_train <- responses[train_idx, ]
  covariates_train <- covariates[train_idx, ]

  n_train <- nrow(responses_train)
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
