library(purrr)
source("cggr_node.R")

#' @param reponses n x d matrix of responses
#' @param covariates n x p matrix of covariates
#' @param lambda_g numeric penalty for mean vector
#' @param lambda_b_seq vector of d penalties for each response
#' @param alpha value between 0 and 1 that determines sparse group lasso penalty
convex_ggr <- function(responses, covariates, lambda_g, lambda_b_seq, alpha,
                       max_iter = 1000,
                       tol = 1e-5) {
  stopifnot(is.matrix(responses),
            is.matrix(covariates),
            nrow(responses) == nrow(covariates),
            lambda_g >= 0,
            all(lambda_b_seq >= 0),
            length(lambda_b_seq) == ncol(responses),
            0 <= alpha && alpha <= 1)

  d <- ncol(responses)
  p <- ncol(covariates)

  # Initialize mean matrix, gamma
  gamma_mx <- matrix(nrow = d, ncol = p)

  # Initialize covariate coef mxs, beta.
  # Includes the population matrix, hence +1
  beta_mxs <-  map(seq_len(p + 1), ~ matrix(nrow = d, ncol = d))

  # Estimated variances
  est_vars <- vector(length = d)

  for (i in seq_len(d)) {
    result <- cggr_node(responses[, i], responses[, -i], covariates,
                        lambda_g, lambda_b_seq[i], alpha, max_iter, tol)[[1]]
    gamma_mx[i, ] <- result$gamma_j
    beta <- result$beta_j

    # Index the columns of each beta matrix to fill. Excludes the diagonal.
    beta_mx_idx <- seq_len(d)[-i]
    for (h in seq_len(p + 1)) {
      # # Set the population mx diagonal entry to the estimated variance, 0 else.
      # diag(beta_mxs[[h]])[i] <- ifelse(h == 1, result$sigma_sq, 0)

      # Set all diagonals to 0. The estimated variance does not really belong in B_0.
      diag(beta_mxs[[h]])[i] <- 0
      bh_idx <- (h - 1) * (d - 1) + seq_len(d - 1)

      # TODO: Is the update to beta correct here?
      # Maybe it is -beta[bh_idx]/result$sigma_sq
      beta_mxs[[h]][i, beta_mx_idx] <- beta[bh_idx]

      est_vars[i] <- result$sigma_sq
    }
  }

  for (h in seq_len(p + 1)) {
    beta_mxs[[h]] <- symmetrize(beta_mxs[[h]])
  }

  return(list(gamma_mx = gamma_mx, beta_mxs = beta_mxs, sigma_sq = est_vars))
}

#' @param covariate p-vector of observation of covariates.
#' @return the mean vector and precision matrix after "undoing" the reparam
est_mvn_params <- function(covariate, result) {
  theta_vec <- result$gamma_mx %*% covariate
  theta_mx <- map2(result$beta_mxs, c(1, as.numeric(covariate)), `*`) |>
                 reduce(`+`)

  diag_prec <- diag(diag(result$beta_mxs[[1]]))
  prec_mx <- - diag_prec %*% theta_mx
  mean_vec <- solve(prec_mx) %*% diag_prec %*% theta_vec

  return(list(mean_vec = mean_vec, prec_mx = prec_mx))
}

#' @return symmetrized version of matrix mx
#' result_ij = result_ji is nonzero iff both mx_ij and mx_ji are nonzero,
#' in which case we choose the smaller value in magnitude.
symmetrize <- function(mx) {
  symm_help <- function(x, y) {
    if (isTRUE(all.equal(x, 0)) & isTRUE(all.equal(y, 0))) {
      return(0)
    }
    ifelse(abs(x) < abs(y), x, y)
  }

  ut <- mx[upper.tri(mx)]
  lt <- t(mx)[upper.tri(mx)]

  symmed <- unlist(map2(ut, lt, symm_help))
  mx[upper.tri(mx)] <- symmed
  for(i in seq_len(nrow(mx))) {
    for(j in seq_len(i - 1)){
      mx[i, j] <- mx[j, i]
    }
  }
  return(mx)
}
