source("convexggr.R")
source("generate_data.R")

set.seed(104)
s <- generate_data(200, 25, 50)
result <- convex_ggr(s$responses[, 1], s$responses[, -1], s$covariates,
                     lambda = c(10, 20), alpha = 0.5, gamma_init = rep(0.25, 50),
                     max_iter = 10)

