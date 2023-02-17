source("test_utils.R")
p <- 25
q <- 100
sparsemean <- FALSE
tbpath <- paste0("data/tb_p", p, "q", q, ".rds")
mgpath <- paste0("data/mg_p", p, "q", q, "sparse", sparsemean, ".rds")
tb <- readRDS(tbpath)
mg <- readRDS(mgpath)
n <- 200

path <- paste0("output/p", p, "q", q, "n", n, ".rds")
cat("Loading results from", path, "\n")
results <- readRDS(path)

avg_perf <- rbind(apply(results$gmmreg, 2, mean),
                  apply(results$gmmreg, 2, sd),
                  apply(results$cggr, 2, mean),
                  apply(results$cggr, 2, sd))
colnames(avg_perf) <- c("TPR", "FPR", "beta_err", "mean_err", "omega_err")
rownames(avg_perf) <- c("gmmreg", "(sd)", "cggr", "(sd)")


for (i in seq_len(6)) {
  cov_lbl <- as.character(i - 1)
  beta_list <- list(GMMReg = results$g_result$bhat[, , i],
                    CGGR = results$c_result$bhat[, , i],
                    Actual = tb[, , i])
  print(beta_viz_list(beta_list, cov_lbl, tileborder = TRUE))
}

print(gamma_viz_list(list(g = results$g_result$ghat, c = results$c_result$ghat, mg)))
