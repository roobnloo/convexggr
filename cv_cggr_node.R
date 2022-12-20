source("utils.R")

cv_cggr_node <- function(node, response, covariates,
                         asparse, regmean, nodereg_strat,
                         nlambda = 100, lambda_factor = 1e-4, nfolds = 5,
                         maxiter = 1000, tol = 1e-5,
                         verbose = FALSE) {
  d <- ncol(response)
  p <- ncol(covariates)
  num <- nrow(response)

  y <- scale(response[, node], scale = FALSE)
  respintx <- cbind(response[, -node], intxmx(response[, -node], covariates))
  lam_max <- max(abs(t(respintx) %*% y))
  lambda <- lam_max * exp(seq(log(1), log(lambda_factor), length = nlambda))

  foldsplit <- split(
    seq_len(num), cut(sample(seq_len(num)), nfolds, labels = FALSE))
  mses <- matrix(nrow = nfolds, ncol = length(lambda))
  for (i in seq_along(foldsplit)) {
    if (verbose)
      print(paste("Begin CV fold", i))
    foldids <- foldsplit[[i]]
    y_train <- y[-foldids]
    response_train <- response[-foldids,]
    covariates_train <- covariates[-foldids,]
    y_test <- y[foldids]
    respintx_test <- respintx[foldids,]
    covariates_test <- covariates[foldids,]

    initbeta <- rep(0, (p + 1) * (d - 1))
    initgamma <- rep(0, p)
    for (lamid in seq_along(lambda)) {
      sgl_node <- nodereg_strat(
        node, response_train, covariates_train,
        lambda[lamid], asparse, regmean,
        initbeta, initgamma, maxit = maxiter, tol = tol)

      resid <- y_test - mean(y_train) - # subtract mean since y is centered
           covariates_test %*% sgl_node$gamma -
           respintx_test %*% sgl_node$beta

      mses[i, lamid] <- sum(resid^2) / (2 * length(resid))
      initgamma <- sgl_node$gamma
      initbeta <- sgl_node$beta
    }
    if (verbose)
      print(paste("Finished CV fold", i))
  }

  cv_mse <- apply(mses, 2, mean)
  return(list(
    cv_mse = cv_mse,
    lambda = lambda,
    lambda_min = lambda[which.min(cv_mse)]
  ))
}
