source("convexggr.R")
source("generate_data.R")

set.seed(104)
d <- 25
p <- 50
num <- 200
lambda = c(1, 5)
s <- generate_data(num, d, p)
result <- convex_ggr(s$responses[, 1], s$responses[, -1], s$covariates,
                     lambda = lambda, alpha = 0.5, gamma_init = rep(0.25, p),
                     max_iter = 50)

