source("convexggr.R")
library(assertthat)

## Test interaction_mx
n <- 5
d <- 3
p <- 3
responses <- matrix(c(rep(1, 5), rep(2, 5)), nrow = n, ncol = 2)
covariates <- matrix(c(rep(3, 5), rep(4, 5), rep(5, 5)), nrow = n, ncol = 3)

test_itx_mx <- interaction_mx(responses, covariates)
expect <- cbind(rep(3, 5), rep(6, 5), rep(4, 5), rep(8, 5), rep(5, 5), rep(10,5))

assert_that(all(test_itx_mx == expect))

## Test sq_error_loss
y_j <- 1:5
gamma_j <- c(0.5, -0.3, 1.4)
b_j0 <- c(2.7, 1.3)
beta_j0 <- seq(1, 5, length = 6)

error_vec <- y_j - covariates %*% gamma_j - responses %*% b_j0 -
             test_itx_mx %*% beta_j0
expected_error <- norm(error_vec, type = "F")^2 / (2 * n)
test_error <- compute_residual(y_j, responses, covariates,
                               gamma_j, c(b_j0, beta_j0)) |>
                sq_error_loss()

assert_that(are_equal(expected_error, test_error))

print("All tests passed!")
