source("gmmreg.R")
source("data_generation.R")

set.seed(100)
p <- 25
q <- 100
tb <- readRDS("data/tb_p25q100.rds")
mg <- readRDS("data/mg_p25q100sparseFALSE.rds")
n <- 200
s <- data_generate(n, tb, mg, reparam = TRUE)

Z <- s$X
U <- s$U
asparse <- seq(0.1, 1, by = 0.3)
nasparse <- length(asparse)
nlambda <- 100

node <- 1

y <- Z[, node]
mx <- cbind(Z[, -node], intxmx(Z[, -node], U))

# There are (q + 1) groups and the size of each group is p-1
grp_idx <- rep(1:(q + 1), each = p - 1)

foldid <- sample(cut(seq_len(n), nfolds, labels = FALSE))
mses <- matrix(nrow = nasparse, ncol = nlambda)

for (i in seq(nasparse)) {
  cv_result <- cv.sparsegl(
    mx, y, grp_idx,
    asparse = asparse[i],
    nlambda = nlambda,
    pf_group = c(0, rep(1, q)),
    foldid = foldid,
    intercept = FALSE,
    standardize = FALSE
  )
  mses[i, ] <- cv_result$cvm
}

# ignore the intercept fitted below with -1
# beta <- as.numeric(coef(cv_result, s = "lambda.min"))[-1]
# rss <- sum((y - predict(cv_result, mx, s = "lambda.min"))^2)
# num_nz <- sum(abs(beta) > 0)
# varhat <- rss / (nrow(Z) - num_nz)
