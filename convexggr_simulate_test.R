library(CVXR)
source("cggr_node.R")
source("generate_data.R")

set.seed(104)
d <- 25
p <- 50
num <- 200
lambda = c(0.1, .3)
s <- generate_data(num, d, p)

# result <- convex_ggr(s$responses, s$covariates,
#                      lambda = lambda, alpha = 0.5,
#                      gamma_init = rep(0.25, p), max_iter = 150)


