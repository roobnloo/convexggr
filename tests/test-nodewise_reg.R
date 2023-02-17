library(Rcpp)
knitr::opts_chunk$set(echo = TRUE)
source("data_generation.R")
source("node_strategy_sparsegl.R")
source("node_strategy_cvxr.R")
source("test_utils.R")
sourceCpp("nodewiseRegression.cpp")


set.seed(104)
p <- 25
q <- 100
tb <- readRDS("data/tb_p25q100.rds")
mg <- readRDS("data/mg_p25q100sparseFALSE.rds")
n <- 200
s <- data_generate(n, tb, mg, reparam = TRUE)

lambda <- 0.19
regmean <- 1.5

result_cpp <- nodewiseRegression(
  s$X[, 1], s$X[, -1], s$U, 0.75,
  regmeanPath = regmean, lambdaPath = lambda,
  maxit = 3000, tol = 1e-8)
cat("done with cpp:", result_cpp$objval, "\n")

result_sgl <- node_strategy_sparsegl(
  s$X[, 1], s$X[, -1], s$U,
  lambda, 0.75, regmean, tol = 1e-8)
result_sgl$objval[length(result_sgl$objval)]
cat("done with sgl:", result_sgl$objval[length(result_sgl$objval)], "\n")

gamma_cvxr <- Variable(q)
beta_cvxr <- Variable((q + 1) * (p - 1))
result_cvxr <- node_strategy_cvxr(
  s$X[, 1], s$X[, -1], s$U, lambda, 0.75,
  regmean, beta_cvxr, gamma_cvxr, tol = 1e-8)
cat("done with cvxr:", result_cvxr$value, "\n")

values <- c(
  result_cvxr$value,
  result_cpp$objval,
  result_sgl$objval[length(result_sgl$objval)]
)

print(values)
print(diff(values))
