library(CVXR)
source("convexggr.R")
source("generate_data.R")

set.seed(104)
d <- 25
p <- 50
num <- 200
lambda = c(0.2, .3)
alpha <- 0.5
s <- generate_data(num, d, p)

my_result <- convex_ggr_component(s$responses[, 1], s$responses[, -1], s$covariates,
                     lambda = lambda, alpha = alpha, max_iter = 1500)

W1 <- interaction_mx(s$responses[, -1], s$covariates)

compute_loss <- function(y, responses, covariates, gamma_j, beta_j) {
  W_j <- interaction_mx(responses, covariates)
  x_gamma_j <- covariates %*% gamma_j # X gamma_j
  y_b_j0 <- responses %*% beta_j[seq_len(d-1)] # Y_-j b_j^0
  W_jbeta_j0 <- W_j %*% beta_j[-seq_len(d-1)] # W_-j beta_j,-0

  r_j <-  y - x_gamma_j - y_b_j0 - W_jbeta_j0
  return(sum_squares(r_j) / (2 * length(r_j)))
}

gamma_j <- Variable(p)
beta_j <- Variable((p + 1) * (d - 1))
loss <- compute_loss(s$responses[, 1], s$responses[, -1], s$covariates,
                     gamma_j, beta_j)

cggr_penalty <- function(gamma_j, beta_j, lambda, alpha){
  group_lasso_term <-
      map(seq_len(p), ~ cvxr_norm(beta_j[.x * (d - 1) + seq_len(d - 1)]), 2) |>
        reduce(`+`)

  lambda[1] * cvxr_norm(gamma_j, p = 1) +
    alpha * lambda[2] * cvxr_norm(beta_j, p = 1) +
    (1 - alpha) * lambda[2] * group_lasso_term
}

obj <- loss + cggr_penalty(gamma_j, beta_j, lambda, alpha)
prob <- Problem(Minimize(obj))
result <- solve(prob)
