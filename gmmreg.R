library(sparsegl)
source("utils.R")

#' @param X n x d matrix of responses
#' @param U n x p matrix of covariates
gmmreg <- function(X, U, asparse, nlambda = 100) {
  stopifnot(is.matrix(X), is.matrix(U),
            nrow(X) == nrow(U),
            asparse > 0, asparse <= 1)
  d <- ncol(X)
  p <- ncol(U)


  # Estimate mean matrix
  ghat_mx <- matrix(nrow = d, ncol = p)
  for (i in seq_len(d)) {
    step1_result <- cv.sparsegl(U, X[, i], seq_len(p),
                                asparse = 1, nfolds = 5,
                                intercept = F, standardize = F)
    # ignore the intercept fitted below with -1
    ghat_mx[i,] <- as.numeric(coef(step1_result, s = "lambda.min"))[-1]
  }

  # Initialize covariate array
  # Includes the population matrix, hence +1
  bhat_tens <-  array(0, dim = c(d, d, p + 1))

  # Estimated variances
  varhat <- vector(length = d)

  Z <- X - U %*% t(ghat_mx)
  for (node in seq_len(d)) {
    y <- scale(Z[, node], scale = F)
    mx <- cbind(Z[, -node], interaction_mx(Z[, -node], U))

    # There are (p + 1) groups and the size of each group is d-1
    grp_idx <- rep(1:(p + 1), each = d-1)

    cv_result <- cv.sparsegl(mx, y, grp_idx,
                             asparse = asparse,
                             nlambda = nlambda,
                             pf_group = c(0, rep(1, p)),
                             nfolds = 5,
                             intercept = F,
                             standardize = F)

    # ignore the intercept fitted below with -1
    beta <- as.numeric(coef(cv_result, s = "lambda.min"))[-1]

    rss <- sum((y - predict(cv_result, mx, s = "lambda.min"))^2)
    num_nz <- sum(abs(beta) > 1e-10)
    varhat[node] <- rss  / (nrow(X) - num_nz)
    bhat_tens[node, -node,] <- -beta / varhat[node]
  }

  bhat_symm <- abind(apply(bhat_tens, 3, symmetrize, simplify = F), along = 3)

  # Returns the estimated precision matrix of the ith observation
  precision <- function(i) {
    omega <- apply(bhat_symm, c(1, 2), \(b) b %*% c(1, U[i, ]))
    diag(omega) <- 1/varhat
    return(omega)
  }

  # Returns the estimated mean vector of the ith observation
  mean <- function(i) {
    mu <- ghat_mx %*% U[i, ]
    return(mu)
  }

  return(list(ghat = ghat_mx,
              bhat = bhat_symm,
              bhat_asym = bhat_tens,
              varhat = varhat,
              precision = precision,
              mean = mean))
}
