library(tidyverse)
library(igraph)

#' @param n number of observations
#' @param p number of responses
#' @param q number of covariates
#' @param qe number of covariates with non-zero effects
#' @param ve Erdos-Renyi edge probability of covariate mxs
#' @param sg number of nonzero mean matrix entries
#' @return a n x p matrix of responses, n x q matrix of covariates,
#' gamma matrix determining the mean, b matrices determining the precision mx,
#' and indices of the covariates with non-zero effects.
generate_data <- function(n, p, q, qe = 5, ve = 0.01, sg = p * q * 0.1) {
  # Generate the true parameters, gamma and beta
  gamma_mx <- generate_mean_mx(p, q, sg)
  b_mxs <- generate_coef_mxs(p, qe, ve)

  # Generate the q covariates
  cov_idx <- seq_len(q)
  cov_nz_idx <- sort(sample(cov_idx, qe)) # indices of non-zero effect covs
  cov_disc_idx <- sort(sample(cov_idx, q/2)) # indices of discrete covs
  covariates <- vector(mode = "list", length = q)
  covariates[cov_disc_idx] <- map(seq_len(q/2), ~ sample(0:1, n, replace = T))
  covariates[-cov_disc_idx] <- map(seq_len(q/2), ~ runif(n))
  covariates <- matrix(unlist(covariates), nrow = n, ncol = q)
  colnames(covariates) <- paste0("u", seq_len(q), sep = "")

  # Generate the p responses
  generate_response <- function(i) {
    mvn_params <- get_mvn_params(i, cov_nz_idx, covariates, gamma_mx, b_mxs)
    response <- rMVNormP(1, mvn_params$mean_vec, mvn_params$prec_mx)
    dim(response) <- NULL
    return(response)
  }
  responses <- map(seq_len(n), generate_response) |>
                unlist() |>
                matrix(nrow = n, ncol = p, byrow = T)
  colnames(responses) <- paste0("x", seq_len(p))

  # Bind p responses and q covariates into the result
  result <- list(responses = responses,
                 covariates = covariates,
                 cov_nz_idx = cov_nz_idx,
                 gamma_mx = gamma_mx,
                 b_mxs = b_mxs)
  return(result)
}

#' @param i index out of n observations
#' @return mean vector and precision matrix for ith response
get_mvn_params <- function(i, cov_nz_idx, covariates, gamma_mx, b_mxs) {

  # Reparametrized mean vector
  little_theta <- gamma_mx %*% as.numeric(covariates[i, ])

  # Reparametrized precision matrix
  big_theta <- map2(b_mxs, c(1, as.numeric(covariates[i, cov_nz_idx])), `*`) |>
               reduce(`+`)

  prec_mx <- -big_theta
  diag(prec_mx) <- 1
  mean_vec <- solve(prec_mx) %*% little_theta

  list(mean_vec = mean_vec, prec_mx = prec_mx)
}

# Generates p x q matrix where sg are non-zero
generate_mean_mx <- function(p, q, sg, val = 0.25) {
  result <- matrix(0, p, q)
  if (sg == 0) {
    return(result)
  }
  idx_sel <- expand.grid(seq_len(p), seq_len(q)) |>
                slice_sample(n = sg)
  result[as.matrix(idx_sel)] <- val
  return(result)
}

# Generate `num` coefs from Unif([-0.5, -0.35] U [0.35, 0.5])
generate_coef <- function(num) {
  neg <- runif(num, min = -0.5, max = -0.35)
  pos <- runif(num, min = 0.35, max = 0.5)
  result_val <- ifelse(rbinom(num, 1, 0.5), neg, pos)
  return(result_val)
}

# Generate one p x p coef matrix for the population-level network
# Underlying graph has scale-free power law node degree
# Fix diagonal entries to -1
# Result is not normalized
generate_pop_coef_mx <- function(p) {
  cov_el <- sample_fitness_pl(p, p, exponent.out = 2.5) |>
              as_edgelist()
  result <- matrix(0, p, p)
  result[cov_el] <- generate_coef(nrow(cov_el))
  diag(result) <- 0
  return(result)
}

# Generate p x p coefficient matrices for each nonzero covariate qe in total
# Underlying graph is ER with edge probability ve.
# Fix diagonals to zero
# Result is not normalized
generate_cov_coef_mxs <- function(p, qe, ve) {
  create_mx <- function(dummy) {
    ntwk_cov <- erdos.renyi.game(p, ve, type = "gnp")
    cov_el <- as_edgelist(ntwk_cov)
    result <- matrix(0, p, p)
    result[cov_el] <- generate_coef(nrow(cov_el))
    diag(result) <- 0

    if (all(result == 0)) {
      return(create_mx(0))
    }
    return(result)
  }
  lapply(seq_len(qe), create_mx)
}

# Generate a list of 1+q coef matrices, including the population-level mx
# Each matrix is normalized to be symmetric and diagonally dominant.
generate_coef_mxs <- function(p, qe, ve) {
  b0_mx <- generate_pop_coef_mx(p)
  b_mxs <- generate_cov_coef_mxs(p, qe, ve)
  result <- c(list(b0_mx), b_mxs)

  # Normalization for diagonal dominance
  for (j in seq_len(p)) {
    normal_j <- sum(unlist(map(result, ~ sum(abs(.x[j, -j])))))
    if (isTRUE(all.equal(normal_j, 0))) {
      next
    }
    for (k in seq_along(result)) {
      result[[k]][j, -j] <- result[[k]][j, -j] / normal_j
    }
  }

  # TODO: clean up this symmetrization
  for (k in seq_along(result)) {
    for (i in seq_len(p)) {
      for (j in seq_len(i-1)) {
        symm <- (result[[k]][i, j] + result[[k]][j, i]) / 2
        result[[k]][i, j] <- -symm
        result[[k]][j, i] <- -symm
      }
    }
  }
  return(result)
}

#' Sample from mvn using precision parametrization
#' @param n number of samples to draw
#' @param mu p-vector mean
#' @param Omega precision matrix (p x p)
#' @return matrix of dimension p x n of samples
rMVNormP <- function(n, mu, Omega){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  U <- chol(Omega)
  X <- backsolve(U, Z) # more efficient and stable than actually inverting
  X <- sweep(X, 1, mu, FUN = `+`)
  return(X)
}
