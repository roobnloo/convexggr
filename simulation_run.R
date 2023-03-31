source("data_generation.R")
source("cggr.R")
source("gmmreg.R")
source("test_utils.R")

# generate_parameters(50, 50, 150, TRUE)

p <- 50
q <- 50
reparam <- TRUE
tbpath <- paste0("data/tb_p", p, "q", q, ".rds")
mgpath <- paste0("data/mg_p", p, "q", q, "sparse", !reparam, ".rds")
tb <- readRDS(tbpath)
mg <- readRDS(mgpath)
n <- 200

generate <- function(n, seed) {
    set.seed(seed)
    data_generate(n, tb, mg, reparam = reparam, verbose = FALSE)
}

failures <- 0
seeds <- 11:50
ntrials <- length(seeds)
s_list <- vector(mode = "list", length = ntrials)
for (i in seq_len(ntrials)) {
  s <- generate(n, seeds[i])
  if (is.null(s)) {
    failures <- failures + 1
    next
  }
  s_list[[i]] <- s
}

cat("Done generating datasets with", failures, "failures\n")

path <- paste0(
  "output/p", p, "q", q, "n", n, "reparam", reparam, "trials", ntrials, ".rds")

gmmreg_result_mx <- matrix(nrow = ntrials, ncol = 5)
cggr_result_mx <- matrix(nrow = ntrials, ncol = 5)

tictoc::tic()
for (i in seq_along(s_list)) {
  g_result <- gmmreg(
    s_list[[i]]$X, s_list[[i]]$U, 0.75,
    verbose = FALSE, parallel = TRUE)
  gmmreg_result_mx[i, ] <- unlist(performance(g_result, s_list[[i]]))

  rm(g_result)
  gc()

  c_result <- cggr(
    s_list[[i]]$X, s_list[[i]]$U, 0.75,
    nregmean = 20,
    verbose = FALSE, parallel = TRUE)
  cggr_result_mx[i, ] <- unlist(performance(c_result, s_list[[i]]))

  rm(c_result)
  gc()

  saveRDS(
    list(gmmreg = gmmreg_result_mx, cggr = cggr_result_mx, seeds = seeds), path)
  cat("Finished round", i, "\n")
}
cat("Finished analysis\n")
tictoc::toc()
cat("Saved results to", path, fill = TRUE)

avg_perf <- rbind(apply(gmmreg_result_mx, 2, mean),
                  apply(gmmreg_result_mx, 2, sd),
                  apply(cggr_result_mx, 2, mean),
                  apply(cggr_result_mx, 2, sd))
colnames(avg_perf) <- c("TPR", "FPR", "beta_err", "mean_err", "omega_err")
rownames(avg_perf) <- c("gmmreg", "(sd)", "cggr", "(sd)")

cat("Average metrics:\n")
print(avg_perf)

# for (i in seq_len(6)) {
#   cov_lbl <- as.character(i - 1)
#   beta_list <- list(GMMReg = g_result$bhat[, , i],
#                     CGGR = c_result$bhat[, , i],
#                     Actual = tb[, , i])
#   print(beta_viz_list(beta_list, cov_lbl, tileborder = TRUE))
# }

# print(gamma_viz_list(list(g = g_result$ghat, c = c_result$ghat, mg)))
