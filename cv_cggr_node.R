library(sparsegl)
source("utils.R")

cv_cggr_node <- function(node, X, U, asparse, regmean,
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
      sgl_node <- cggr_node(node, X_train, U_train, lambda[lamid], asparse,
                            initbeta, initgamma,
                            regmean = regmean, maxiter = maxiter, tol = tol)
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

cggr_node <- function(node, X, U, lambda, asparse, regmean,
                      initbeta = NULL, initgamma = NULL,
                      maxiter = 1000, tol = 1e-5) {
  d <- ncol(X)
  p <- ncol(U)
  n <- nrow(X)

  stopifnot(1 <= node && node <= d,
            nrow(U) == n,
            asparse >= 0, length(asparse) == 1,
            lambda >= 0, length(lambda) == 1)

  y <- scale(X[, node], scale = F)
  W <- cbind(X[, -node], interaction_mx(X[, -node], U))

  if (is.null(initbeta)) {
    initbeta <- rep(0, (p+1)*(d-1))
  }
  if (is.null(initgamma)) {
    initgamma <- rep(0, p)
  }
  beta <- initbeta
  gamma <- initgamma

  # There are (p + 1) groups and the size of each group is d-1
  grp_idx <- rep(1:(p + 1), each = d-1)

  # residual vector
  r <- y - U %*% gamma - W %*% beta
  objval <- vector(length = maxiter) # keep track of objective values

  for (i in seq_len(maxiter)) {
    objval[i] <- sum(r^2) / (2 * n) + regmean * sum(gamma^2) +
                   asparse * lambda * sum(abs(beta[seq_len(d-1)])) +
                   (1 - asparse) * lambda * sum(beta[-seq_len(d-1)]^2)
    if (i > 1 && abs(objval[i] - objval[i-1]) < tol) {
      objval <- objval[1:i]
      break
    }

    # sparsegl step
    rbeta <- r + W %*% beta # partial residual for beta
    sgl <- sparsegl(W, rbeta, grp_idx, lambda = lambda, asparse = asparse,
                    intercept = F, standardize = F, pf_group = c(0, rep(1, p)))
    beta <- as.numeric(coef(sgl, s = lambda)[-1])
    r <- rbeta - W %*% beta

    # ridge step
    rgamma <- r + U %*% gamma
    gamma <- solve(t(U) %*% U + diag(regmean, p, p)) %*% t(U) %*% rgamma
    r <- rgamma - U %*% gamma
  }

  return(list(beta = beta,
              gamma = gamma,
              objval = objval,
              resid = r))
}
