test <- function() {
  source("data_generation.R")
  source("cggr.R")
  source("test_utils.R")

  set.seed(101)
  n <- 100
  p <- 5
  q <- 10
  dim <- (p - 1) * (q + 1)
  # load("data/tB.RData")
  tB <- generate_tb(p, q)
  s <- data_generate(n, p, q, tB)

  cggr(s$X, s$U, asparse = 0.75)
}

result <- test()
