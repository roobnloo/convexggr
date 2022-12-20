source("utils.R")
library(sparsegl)

node_strategy_sparsegl <- function(node, X, U, lambda, asparse, regmean,
                                   initbeta = NULL, initgamma = NULL,
                                   maxit = 1000, tol = 1e-5) {
  d <- ncol(X)
  p <- ncol(U)
  n <- nrow(X)

  stopifnot(1 <= node && node <= d,
            nrow(U) == n,
            asparse >= 0, length(asparse) == 1,
            lambda >= 0, length(lambda) == 1)

  y <- scale(X[, node], scale = F)
  W <- cbind(X[, -node], intxmx(X[, -node], U))

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
  objval <- vector(length = maxit) # keep track of objective values

  for (i in seq_len(maxit)) {
    objval[i] <- sum(r^2) / (2 * n) + regmean * sum(gamma^2) +
                   asparse * lambda * sum(abs(beta)) +
                   (1 - asparse) * lambda *
                    sum(sapply(seq_len(p),
                             \(x) sqrt(sum((beta[x * (d-1) + seq_len(d-1)])^2)))
                    )
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
    gamma <- solve(t(U) %*% U + diag(regmean, p, p), t(U) %*% rgamma)
    r <- rgamma - U %*% gamma
  }

  return(list(beta = beta,
              gamma = gamma,
              objval = objval,
              resid = r))
}
