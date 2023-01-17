source("utils.R")

cv_cggr_node <- function(node, responses, covariates,
                         lambdas, asparse, regmean,
                         maxit, tol, nfolds,
                         verbose = FALSE, parallel = FALSE) {
  p <- ncol(responses)
  q <- ncol(covariates)
  n <- nrow(responses)
  nlambda <- length(lambdas)

  intx <- cbind(responses[, -node], intxmx(responses[, -node], covariates))
  foldsplit <- split(
    seq_len(n), cut(sample(seq_len(n)), nfolds, labels = FALSE))
  mses <- matrix(nrow = nfolds, ncol = nlambda)
  for (i in seq_len(nfolds)) {
    testids <- foldsplit[[i]]
    y_train <- responses[-testids, node]
    responses_train <- responses[-testids, -node]
    covariates_train <- covariates[-testids, ]

    ntest <- length(testids)
    y_test <- responses[testids, node]
    y_test <- matrix(rep(y_test, times = nlambda), nrow = ntest, ncol = nlambda)
    intx_test <- intx[testids, ]
    covariates_test <- covariates[testids, ]

    nodereg <- nodewiseRegression(
      y_train, responses_train, covariates_train, asparse, regmean,
      lambdas = lambdas, maxit = maxit, tol = tol)

    resid_test <- y_test - covariates_test %*% nodereg["gamma"][[1]] -
                  intx_test %*% nodereg["beta"][[1]]
    mses[i, ] <- apply(resid_test, 2, \(x) sum(x^2) / (ntest - 1))
    # if (verbose)
    #   print(paste("Finished CV fold", i))
  }

  cv_mse <- apply(mses, 2, mean)
  return(list(
    cv_mse = cv_mse,
    lambda_min = lambdas[which.min(cv_mse)]
  ))
}
