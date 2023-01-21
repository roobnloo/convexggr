library(Rcpp)
knitr::opts_chunk$set(echo = TRUE)
source("data_generation.R")
source("cggr.R")
# source("node_strategy_sparsegl.R")
# source("node_strategy_cvxr.R")
source("test_utils.R")
# sourceCpp("nodewiseRegression.cpp")


set.seed(104)
n <- 200
p <- 10
q <- 20
mg <- generate_mg(p, q)
tb <- generate_tb(p, q, 0.02)
# tb <- readRDS("data/tB.rds")
s <- data_generate(n, p, q, tb, mg, reparam = FALSE)

# nodewiseRegression(
#   s$X[, 1], s$X[, -1], s$U, 0.20, 0.01,
#   maxit = 1000, tol = 1e-6, verbose = TRUE)
# print("goodbye")

tictoc::tic()
result_cpp <- cggr(
  s$X, s$U, 0.20, 0.01, tol = 1e-5, verbose = TRUE, parallel = TRUE)
tictoc::toc()
result_cpp$lambda[cbind(result_cpp$cv_lambda, 1:p)]

# tictoc::tic()
# result_sgl <- cggr(
#   s$X, s$U, 0.20, node_strategy_sparsegl, 0.01, maxit = 1000, tol = 1e-5,
#   verbose = TRUE)
# tictoc::toc()
# # print(result_cpp$objval)
# result_sgl$lambda[cbind(result_sgl$cv_lambda, 1:p)]
# print(gamma_viz_compare(result_sgl$ghat, mG))
# print(beta_viz_compare(result_sgl$bhat_asym[, , 1], tb[, , 1],
#       "population", tileborder = TRUE))

# tictoc::tic()
# result_sgl <- node_strategy_sparsegl(
#   s$X[, 1], s$X[, -1], s$U,
#   0.03, 0.75, 0.01, initbeta, initgamma, tol = 1e-8)
# result_sgl$objval[length(result_sgl$objval)]
# tictoc::toc()

# gamma_cvxr <- Variable(q)
# beta_cvxr <- Variable((q + 1) * (p - 1))
# result_cvxr <- node_strategy_cvxr(
#   s$X[, 1], s$X[, -1], s$U,
#   result_cpp$lambdas[9], 0.75, 0.01, beta_cvxr, gamma_cvxr, maxit = 5000, tol = 1e-8)
# values <- c(
#   result_cvxr$value,
#   result_cpp$objval[length(result_cpp$objval)],
#   result_sgl$objval[length(result_sgl$objval)]
# )
# values
