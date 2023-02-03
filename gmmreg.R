library(abind)
library(sparsegl)
source("utils.R")

#' @param X n x d matrix of responses
#' @param U n x p matrix of covariates
gmmreg <- function(X, U, asparse, nlambda = 100, verbose = FALSE, parallel = TRUE) {
  stopifnot(is.matrix(X), is.matrix(U),
            nrow(X) == nrow(U),
            asparse > 0, asparse <= 1)
  p <- ncol(X)
  q <- ncol(U)


  # Estimate mean matrix
  ghat_mx <- matrix(nrow = p, ncol = q)

  nodewise_gamma <- function(node) {
    if (verbose)
      print("Begin step 1 for node", node, quote = FALSE)
    result <- cv.sparsegl(U, X[, node], seq_len(q),
                                asparse = 1, nfolds = 5,
                                intercept = FALSE, standardize = FALSE)
    # print(paste("lambda min for node", node, "is", result$lambda.min))
    if (verbose)
      print("Finished step 1 for node", node, quote = FALSE)
    # ignore the intercept fitted below with -1
    gamma <- as.numeric(coef(result, s = "lambda.min"))[-1]
    return(gamma)
  }

  if (parallel) {
    step1_result <- parallel::mclapply(seq_len(p), nodewise_gamma, mc.cores = 13L);
  } else {
    step1_result <- lapply(seq_len(p), nodewise_gamma);
  }

  for (node in seq_len(p)) {
    ghat_mx[node, ] <- step1_result[[node]]
  }

  # Initialize covariate array
  # Includes the population matrix, hence +1
  bhat_tens <-  array(0, dim = c(p, p, q + 1))

  # Estimated variances
  varhat <- vector(length = p)

  Z <- X - U %*% t(ghat_mx)

  nodewise_beta <- function(node)
  {
    y <- Z[, node]
    mx <- cbind(Z[, -node], intxmx(Z[, -node], U))

    # There are (q + 1) groups and the size of each group is p-1
    grp_idx <- rep(1:(q + 1), each = p - 1)

    if (verbose)
      print("Begin regression for node", node, quote = FALSE)
    cv_result <- cv.sparsegl(mx, y, grp_idx,
      asparse = asparse,
      nlambda = nlambda,
      pf_group = c(0, rep(1, q)),
      nfolds = 5,
      intercept = FALSE,
      standardize = FALSE
    )
    if (verbose)
      print("Finished regression for node", node, quote = FALSE)

    # ignore the intercept fitted below with -1
    beta <- as.numeric(coef(cv_result, s = "lambda.min"))[-1]
    rss <- sum((y - predict(cv_result, mx, s = "lambda.min"))^2)
    num_nz <- sum(abs(beta) > 0)
    return(list(
      beta = beta,
      varhat = rss  / (nrow(X) - num_nz)
    ))
  }

  if (parallel) {
    result <- parallel::mclapply(seq_len(p), nodewise_beta, mc.cores = 13L)
  } else {
    result <- lapply(seq_len(p), nodewise_beta)
  }
  for (node in seq_len(p)) {
    varhat[node] <- result[[node]]$varhat
    bhat_tens[node, -node,] <- result[[node]]$beta / varhat[node]
  }

  bhat_symm <- abind(apply(bhat_tens, 3, symmetrize, simplify = F), along = 3)

  # Returns the estimated precision matrix of the ith observation
  precision <- function(i) {
    omega <- -apply(bhat_symm, c(1, 2), \(b) b %*% c(1, U[i, ]))
    diag(omega) <- 1 / varhat
    return(omega)
  }

  # Returns the estimated mean vector of the ith observation
  mean <- function(i) {
    mu <- ghat_mx %*% U[i, ]
    return(mu)
  }

  list(
    ghat = ghat_mx,
    bhat = bhat_symm,
    bhat_asym = bhat_tens,
    varhat = varhat,
    precision = precision,
    mean = mean
  )
}
