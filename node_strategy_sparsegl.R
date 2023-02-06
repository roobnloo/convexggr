source("utils.R")
library(sparsegl)

node_strategy_sparsegl <- function(y, responses, covariates,
                                   lambda, asparse, regmean,
                                   initbeta = NULL, initgamma = NULL,
                                   maxit = 1000, tol = 1e-5) {
  p <- ncol(responses) + 1
  q <- ncol(covariates)
  n <- nrow(responses)

  W <- cbind(responses, intxmx(responses, covariates))

  stopifnot(nrow(covariates) == n,
            asparse >= 0, length(asparse) == 1,
            lambda >= 0, length(lambda) == 1)

  if (is.null(initbeta)) {
    initbeta <- rep(0, (q + 1) * (p - 1))
  }
  if (is.null(initgamma)) {
    initgamma <- rep(0, q)
  }
  beta <- initbeta
  gamma <- initgamma

  # There are (q + 1) groups and the size of each group is p-1
  grp_idx <- rep(1:(q + 1), each = p-1)

  # residual vector
  r <- y - covariates %*% gamma - W %*% beta
  objval <- vector(length = maxit) # keep track of objective values

  for (i in seq_len(maxit)) {
    objval[i] <- sum(r^2) / (2 * n) + regmean * sum(gamma^2) +
                   asparse * lambda * sum(abs(beta)) +
                   (1 - asparse) * lambda *
                    sum(sapply(seq_len(q),
                             \(x) sqrt(sum((beta[x * (p-1) + seq_len(p-1)])^2)))
                    )
    if (i > 1 && abs(objval[i] - objval[i-1]) < tol) {
      objval <- objval[1:i]
      break
    }

    # sparsegl step
    rbeta <- r + W %*% beta # partial residual for beta
    sgl <- sparsegl(W, rbeta, grp_idx, lambda = lambda, asparse = asparse,
                    intercept = F, standardize = F, pf_group = c(0, rep(1, q)))
    beta <- as.numeric(coef(sgl, s = lambda)[-1])
    r <- rbeta - W %*% beta

    # ridge step
    rgamma <- r + covariates %*% gamma
    gamma <- solve(
      t(covariates) %*% covariates / n + diag(regmean * 2, q, q),
      t(covariates) %*% rgamma) / n
    r <- rgamma - covariates %*% gamma
  }

  return(list(beta = beta,
              gamma = gamma,
              objval = objval,
              resid = r))
}
