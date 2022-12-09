library(sparsegl)
source("utils.R")

cv_cggr_node <- function(node, X, U, asparse, regmean, nodereg_strat,
                         nlambda = 100, lambda_factor = 1e-4, nfolds = 5,
                         maxiter = 1000, tol = 1e-5) {
  d <- ncol(X)
  p <- ncol(U)
  num <- nrow(X)

  y <- scale(X[, node], scale = F)
  W <- cbind(X[, -node], interaction_mx(X[, -node], U))
  lam_max <- max(abs(t(W) %*% y))
  lambda <- lam_max * exp(seq(log(1), log(lambda_factor), length = nlambda))

  foldsplit <- split(seq_len(num), cut(sample(seq_len(num)), nfolds, labels = F))
  mses <- matrix(nrow = nfolds, ncol = length(lambda))
  for(i in seq_along(foldsplit)) {
    foldids <- foldsplit[[i]]
    y_train <- y[-foldids]
    X_train <- X[-foldids,]
    U_train <- U[-foldids,]
    y_test <- y[foldids]
    W_test <- W[foldids,]
    U_test <- U[foldids,]

    initbeta <- rep(0, (p+1)*(d-1))
    initgamma <- rep(0, p)
    for(lamid in seq_along(lambda)) {
      sgl_node <- nodereg_strat(node, X_train, U_train, lambda[lamid],
                                asparse, initbeta, initgamma,
                                regmean = regmean, maxit = maxiter, tol = tol)
      r <- y_test - mean(y_train) - # subtract mean since y is centered
           U_test %*% sgl_node$gamma -
           W_test %*% sgl_node$beta
      mses[i, lamid] <- sum(r^2) / (2*length(r))
      initgamma <- sgl_node$gamma
      initbeta <- sgl_node$beta
    }
  }

  cv_mse <- apply(mses, 2, mean)
  return(list(
    cv_mse = cv_mse,
    lambda = lambda,
    lambda_min = lambda[which.min(cv_mse)]
  ))
}
