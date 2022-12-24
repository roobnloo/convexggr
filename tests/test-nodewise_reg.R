library(Rcpp)
knitr::opts_chunk$set(echo = TRUE)
source("data_generation.R")
source("cggr.R")
source("node_strategy_sparsegl.R")
source("node_strategy_custom.R")
source("node_strategy_cvxr.R")
source("test_utils.R")
sourceCpp("nodewise_regression.cpp")

set.seed(102)
n <- 200
p <- 50
q <- 50
mg <- generate_mg(p, q)
tb <- generate_tb(p, q, 0.02)
# tb <- readRDS("data/tB.rds")
s <- data_generate(n, p, q, tb, mg, reparam = FALSE)
initbeta <- rep(0, (p - 1) * (q + 1))
initgamma <- rep(0, q)

tictoc::tic()
result_cpp <- nodewiseRegression(
    s$X[, 1], s$X[, -1], s$U,
    0.03, 0.75, 0.01, initbeta, initgamma, tol = 1e-10)
tictoc::toc()

result_sgl <- node_strategy_sparsegl(
  s$X[, 1], s$X[, -1], s$U,
  0.03, 0.75, 0.01, initbeta, initgamma, tol = 1e-10)
result_sgl$objval[length(result_sgl$objval)]

gamma_cvxr <- Variable(q)
beta_cvxr <- Variable((q + 1) * (p - 1))
result_cvxr <- node_strategy_cvxr(
  s$X[, 1], s$X[, -1], s$U,
  0.03, 0.75, 0.01, beta_cvxr, gamma_cvxr, maxit = 5000, tol = 1e-12)
values <- c(
  result_cvxr$value,
  result_cpp$objval[length(result_cpp$objval)],
  result_sgl$objval[length(result_sgl$objval)]
)
values
