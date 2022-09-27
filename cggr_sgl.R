library(sparsegl)
source("utils.R")

#' @param reponses n x d matrix of responses
#' @param covariates n x p matrix of covariates
#' @param lambda_g numeric penalty for mean vector
#' @param lambda_b_seq vector of d penalties for each response
#' @param alpha_seq values between 0 and 1 that determines sparse group lasso penalty
cggr_sgl <- function(responses, covariates, alpha_seq, reg_gamma = 0.01,
                     reparametrize = T) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            all(alpha_seq >= 0, alpha_seq <= 1))

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
    mx <- cbind(responses[, -node],
                interaction_mx(responses[, -node], covariates))

    # There are (p + 1) groups and the size of each group is d-1
    grp_idx <- rep(1:(p + 1), each = d-1)
    cv_result <- cv.sparsegl(mx, responses[, node], grp_idx,
                             asparse = alpha_seq,
                             # pf_group = c(0, rep(sqrt(d-1), p)),
                             pf_group = c(0, rep(1, p)),
                             intercept = F)

    # gamma_mx[node, ] <- result$gamma_j
    # ignore the intercept fitted below with -1
    beta <- as.numeric(coef(cv_result, s = "lambda.min"))[-1]

    rss <- sum((responses[, node] - predict(cv_result, mx,s = "lambda.min"))^2)
    num_nz <- sum(abs(beta) > 1e-10) # add gamma to this later
    est_vars[node] <- rss  / (nrow(responses) - num_nz)

    # Index the columns of each beta matrix to fill.
    # Diagonals are set to zero.
    beta_mx_idx <- seq_len(d)[-node]
    for (h in seq_len(p + 1)) {
      diag(beta_mxs[[h]])[node] <- 0
      bh_idx <- (h - 1) * (d - 1) + seq_len(d - 1)

      if (reparametrize) {
        beta_mxs[[h]][node, beta_mx_idx] <- beta[bh_idx]
      } else {
        beta_mxs[[h]][node, beta_mx_idx] <- -beta[bh_idx] / est_vars[node]
      }
    }
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
