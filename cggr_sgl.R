library(sparsegl)
source("utils.R")

#' @param reponses n x d matrix of responses
#' @param covariates n x p matrix of covariates
#' @param lambda_g numeric penalty for mean vector
#' @param lambda_b_seq vector of d penalties for each response
#' @param alpha_seq values between 0 and 1 that determines sparse group lasso penalty
cggr_sgl <- function(responses, covariates, alpha_seq, regmean = 0.01,
                     reparametrize = T, lambda_seq = NULL, tune = T) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            all(alpha_seq >= 0, alpha_seq <= 1),
            tune || !is.null(lambda_seq))

  d <- ncol(responses)
  p <- ncol(covariates)

  # Initialize mean matrix, gamma
  gamma_mx <- matrix(nrow = d, ncol = p)

  # Initialize covariate coef mxs, beta.
  # Includes the population matrix, hence +1
  beta_mxs <-  map(seq_len(p + 1), ~ matrix(nrow = d, ncol = d))

  # Estimated variances
  est_vars <- vector(length = d)

  for (node in seq_len(d)) {
    result <- cggr_sgl_node(node, responses, covariates, lambda = 0.03, alpha = 0.75)
    rss <- sum(result$resid^2)
    num_nz <- sum(abs(result$beta) > 1e-10) + sum(abs(result$gamma) > 1e-10)
    est_vars[node] <- rss  / (nrow(responses) - num_nz)

    # Index the columns of each beta matrix to fill.
    # Diagonals are set to zero.
    beta_mx_idx <- seq_len(d)[-node]
    for (h in seq_len(p + 1)) {
      diag(beta_mxs[[h]])[node] <- 0
      bh_idx <- (h - 1) * (d - 1) + seq_len(d - 1)

      if (reparametrize) {
        beta_mxs[[h]][node, beta_mx_idx] <- result$beta[bh_idx]
      } else {
        beta_mxs[[h]][node, beta_mx_idx] <- -result$beta[bh_idx] / est_vars[node]
      }
    }

    gamma_mx[node,] <- result$gamma
  }

  symmed_beta_mxs <-  map(seq_len(p + 1), ~ matrix(nrow = d, ncol = d))
  for (h in seq_len(p + 1)) {
    symmed_beta_mxs[[h]] <- symmetrize(beta_mxs[[h]])
  }

  return(list(gamma_mx = gamma_mx,
              beta_mxs = symmed_beta_mxs,
              asym_beta_mx = beta_mxs,
              sigma_sq = est_vars))
}

cv.cggr_sgl_node <- function(node, resp, covar, nlambda, lambda_factor = 1e-4,
                             asparse,
                             initbeta = NULL, initgamma = NULL,
                             nfold = 5,
                             regmean = 0.01,
                             maxiter = 1000, tol = 1e-5) {
  d <- ncol(resp)
  p <- ncol(covar)
  num <- nrow(resp)
  beta <- initbeta
  gamma <- initgamma
  if (is.null(beta)) {
    beta <- rep(0, (p+1)*(d-1))
  }
  if (is.null(gamma)) {
    gamma <- rep(0, p)
  }

  y <- scale(resp[, node], scale = F)
  W <- cbind(resp[, -node], interaction_mx(resp[, -node], covar))

  lam_max <- max(abs(t(W) %*% y))
  lambda_seq <- lam_max * exp(seq(log(1), log(lambda_factor)))

  foldids <- split(seq_len(num), cut(seq_along(num), nfold, labels = F))
}

cggr_sgl_node <- function(node, resp, covar, lambda, alpha,
                          initbeta = NULL, initgamma = NULL,
                          regmean = 0.01,
                          maxiter = 1000, tol = 1e-5) {
  d <- ncol(resp)
  p <- ncol(covar)
  num <- nrow(resp)
  beta <- initbeta
  gamma <- initgamma
  if (is.null(beta)) {
    beta <- rep(0, (p+1)*(d-1))
  }
  if (is.null(gamma)) {
    gamma <- rep(0, p)
  }

  y <- scale(resp[, node], scale = F)
  W <- cbind(resp[, -node], interaction_mx(resp[, -node], covar))

  # There are (p + 1) groups and the size of each group is d-1
  grp_idx <- rep(1:(p + 1), each = d-1)

  # residual vector
  r <- y - covar %*% gamma - W %*% beta
  objval <- vector(length = 1000) # keep track fo objective values

  for (i in seq_len(maxiter)) {
    objval[i] <- sum(r^2) / (2 * num) + regmean * sum(gamma^2) +
                   alpha * lambda * sum(abs(beta[seq_len(d-1)])) +
                   (1 - alpha) * lambda * sum(beta[-seq_len(d-1)]^2)
    if (i > 1 && abs(objval[i] - objval[i-1]) < tol) {
      objval <- objval[1:i]
      break
    }

    # sparsegl step
    rbeta <- r + W %*% beta # partial residual for beta
    sgl <- sparsegl(W, rbeta, grp_idx, lambda = lambda, asparse = alpha,
                    intercept = F, standardize = F, pf_group = c(0, rep(1, p)))
    beta <- as.numeric(coef(sgl, s = lambda)[-1])
    r <- rbeta - W %*% beta

    # ridge step
    rgamma <- r + covar %*% gamma
    gamma <- solve(t(covar) %*% covar + diag(regmean, p, p)) %*% t(covar) %*% rgamma
    r <- rgamma - covar %*% gamma
  }

  return(list(beta = beta,
              gamma = gamma,
              objval = objval,
              resid = r))
}
