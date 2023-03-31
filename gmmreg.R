library(abind)
library(sparsegl)
source("utils.R")

#' @param responses n x p matrix of responses
#' @param covariates n x q matrix of covariates
gmmreg <- function(
  responses, covariates, asparse = seq(0.1, 1, by = 0.1), nlambda = 100,
  nfolds = 5, verbose = FALSE, parallel = TRUE) {
  stopifnot(is.matrix(responses), is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            all(asparse >= 0), all(asparse <= 1))
  p <- ncol(responses)
  q <- ncol(covariates)
  n <- nrow(responses)
  nasparse <- length(asparse)

  # Estimate mean matrix
  ghat_mx <- matrix(nrow = p, ncol = q)

  nodewise_gamma <- function(node) {
    result <- cv.sparsegl(covariates, responses[, node], seq_len(q),
                                asparse = 1, nfolds = nfolds,
                                intercept = FALSE, standardize = FALSE)
    # print(paste("lambda min for node", node, "is", result$lambda.min))
    if (verbose)
      cat("Finished step 1 for node", node, fill = TRUE)
    # ignore the intercept fitted below with -1
    gamma <- as.numeric(coef(result, s = "lambda.min"))[-1]
    return(gamma)
  }

  if (parallel) {
    step1_result <- parallel::mclapply(seq_len(p), nodewise_gamma, mc.cores = 10L)
  } else {
    step1_result <- lapply(seq_len(p), nodewise_gamma)
  }

  for (node in seq_len(p)) {
    ghat_mx[node, ] <- step1_result[[node]]
  }

  rm(step1_result)
  gc()

  # Initialize covariate array
  # Includes the population matrix, hence +1
  bhat_tens <-  array(0, dim = c(p, p, q + 1))

  # Estimated variances
  varhat <- vector(length = p)

  Z <- responses - covariates %*% t(ghat_mx)

  nodewise_beta <- function(node) {
    y <- Z[, node]
    mx <- cbind(Z[, -node], intxmx(Z[, -node], covariates))

    # There are (q + 1) groups and the size of each group is p-1
    grp_idx <- rep(1:(q + 1), each = p - 1)

    # if (verbose)
    #   cat("Begin regression for node", node, fill = TRUE)

    foldid <- sample(cut(seq_len(n), nfolds, labels = FALSE))
    mses <- matrix(nrow = nlambda, ncol = nasparse)
    betas <- array(dim = c((p - 1) * (q + 1), nlambda, nasparse))
    varhat <- matrix(nrow = nlambda, ncol = nasparse)

    for (i in seq_len(nasparse)) {
      cv_result <- cv.sparsegl(
        mx, y, grp_idx,
        asparse = asparse[i],
        nlambda = nlambda,
        pf_group = c(0, rep(1, q)),
        foldid = foldid,
        intercept = FALSE,
        standardize = FALSE
      )
      mses[, i] <- cv_result$cvm
      # -1 ignores the intercept
      betas[, , i] <- as.numeric(coef(cv_result, s = cv_result$lambda)[-1, ])
      nnzero_pop <- apply(betas[seq_len(p - 1), , i], 2, \(b) sum(abs(b) > 0))
      # Adjusts for double-counting if population group is active
      nnzero_pop[nnzero_pop == 0] <- 1
      nnzero_pop <- nnzero_pop - 1
      rss <- apply((y - predict(cv_result, mx, s = cv_result$lambda))^2, 2, sum)
      varhat[, i] <- rss / (n - (cv_result$active_grps + nnzero_pop))
    }
    cvind <- arrayInd(which.min(mses), dim(mses))
    if (verbose) {
      cat("Finished regression for node", node, fill = TRUE)
      cat("CV indices:", sprintf("(%d, %d)", cvind[1], cvind[2]), fill = TRUE)
    }

    if (length(varhat[cvind[1], cvind[2]]) == 0) {
      cat("CV indices:", sprintf("(%d, %d)", cvind[1], cvind[2]), fill = TRUE)
      print(varhat)
    }

    return(list(
      beta = betas[, cvind[1], cvind[2]],
      varhat = varhat[cvind[1], cvind[2]],
      cvind = cvind
    ))
  }

  if (parallel) {
    result <- parallel::mclapply(seq_len(p), nodewise_beta, mc.cores = 10L)
  } else {
    result <- lapply(seq_len(p), nodewise_beta)
  }
  for (node in seq_len(p)) {
    varhat[node] <- result[[node]]$varhat
    bhat_tens[node, -node, ] <- result[[node]]$beta / varhat[node]
  }

  bhat_symm <- abind(
    apply(bhat_tens, 3, symmetrize, simplify = FALSE), along = 3)

  # Returns the estimated precision matrix of the ith observation
  precision <- function(i) {
    omega <- -apply(bhat_symm, c(1, 2), \(b) b %*% c(1, covariates[i, ]))
    diag(omega) <- 1 / varhat
    return(omega)
  }

  # Returns the estimated mean vector of the ith observation
  mean <- function(i) {
    mu <- ghat_mx %*% covariates[i, ]
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
