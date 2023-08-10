library(MASS)
library(igraph)

generate_parameters <- function(p, q, seed, sg = 125, sparse_mean = TRUE) {
    set.seed(seed)
    tb <- generate_tb(p, q)
    if (sparse_mean) {
        mg <- generate_mg(p, q, sg)
    } else {
        mg <- matrix(runif(p * q, -0.25, 0.25), p, q)
    }
    tbpath = paste0("data/tb_p", p, "q", q, ".rds")
    mgpath = paste0("data/mg_p", p, "q", q, "sparse", sparse_mean, ".rds")
    saveRDS(tb, tbpath)
    saveRDS(mg, mgpath)
    cat("Done saving.")
}

data_generate <- function(n, tB, mG, reparam = FALSE, verbose = TRUE) {
  pd_test <- 0
  try <- 0
  while (pd_test == 0) {
    output <- X_simulate(n, tB, mG, reparam)
    pd_test <- min(output$pd_check)
    #  if this value is greater than 0, it means all covariances are PD
    if (verbose)
      print(pd_test)

    if (try == 10) {
      print(paste("Failed after", try, "tries"))
      return()
    }
    try <- try + 1
  }
  if (verbose)
    print(paste("Completed after", try, "tries."))
  return(output)
}

X_simulate <- function(n, tB, mG, reparam) {
  p <- nrow(mG)
  q <- dim(tB)[3] - 1
  U <- matrix(sample(c(0, 1), n * q, replace = TRUE), n, q)
  U2 <- matrix(runif(n * q), n, q)
  cind <- sample(1:q, q / 2, replace = FALSE)
  U[, cind] <- U2[, cind] # create covariates that are continuous and binary
  U <- apply(U, 2, scale)
  iU <- cbind(rep(1, n), U)
  X <- matrix(0, n, p)
  mumx <- matrix(0, n, p) # mean vector for each observation

  for (i in 1:n) {
    omega <- -apply(tB, c(1, 2), \(b) b %*% iU[i, ])
    diag(omega) <- 1
    sigma <- solve(omega)
    if (!all(eigen(omega)$values > 0)) { # break if omega(u) is not PD
      break
    }
    mu <- mG %*% U[i,]
    if (reparam) {
      mu <- sigma %*% mu
    }
    mumx[i,] <- mu
    X[i, ] <- mvrnorm(1, mu, sigma)
  }
  return(list(X = X, U= U,
              mG = mG, tB = tB,
              mumx = mumx,
              pd_check = apply(abs(X), 1, sum)))
}

generate_mg <- function(p, q, sg = 125) {
  G <- matrix(0, nrow = p, ncol = q)
  G[sample(seq_along(G), sg)] <- 0.25
  return(G)
}

generate_tb <- function(p, q, ve = 0.01, nz_cov = 5) {
  tB <- array(0, dim = c(p, p, q + 1)) #+1 for intercept
  l <- 0.35
  u <- 0.5

  # degs<-c(2,2,1,1,1,1,rep(0,p-6))
  # g <- sample_degseq(degs, method="simple.no.multiple") #generate with given degrees

  g1 <- sample_pa(p, power = 1.1, directed = FALSE) # scale-free network
  A <- get.adjacency(g1, sparse = F)
  rind <- sample(1:p, p, replace = FALSE)
  A <- A[rind, rind]
  tb <- matrix(0, p, p)
  tb[lower.tri(A) & A > 0] <- sample(c(runif(sum(A), -u, -l),
                                       runif(sum(A), l, u)),
                                     sum(A) / 2, replace = FALSE)
  tB[, , 1] <- tb + t(tb)

  for (j in seq(2, nz_cov + 1)) {
    g2 <- sample_gnp(p, ve, directed = FALSE, loops = FALSE) # random network
    A <- get.adjacency(g2, sparse = F)
    rind <- sample(1:p, p, replace = FALSE)
    A <- A[rind, rind]
    tb <- matrix(0, p, p)
    tb[lower.tri(A) & A > 0] <- sample(c(runif(sum(A), -u, -l),
                                         runif(sum(A), l, u)),
                                       sum(A) / 2, replace = F)
    tB[, , j] <- tb + t(tb)
  }

  tB_temp <- array(0, dim = c(p, p, q + 1))
  for (j in 1:p) {
    tB_temp[, j, ] <- tB[, j, ] / (sum(abs(tB[, j, ])) * 1.5) # ensure diagonal dominance
  }

  for (j in 1:(q + 1)) {
    tB[, , j] <- (tB_temp[, , j] + t(tB_temp[, , j])) / 2
  }

  return(tB)
}
