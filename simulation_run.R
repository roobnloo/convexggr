source("data_generation.R")
source("cggr.R")
source("gmmreg.R")
source("test_utils.R")

# generate_parameters(50, 50, 150, TRUE)

p <- 50
q <- 100
reparam <- TRUE
sparse <- TRUE
sg <- 500
tbpath <- paste0("data/tb_p", p, "q", q, ".rds")
mgpath <- paste0("data/mg_p", p, "q", q, "sparse", sg, ".rds")
tb <- readRDS(tbpath)
mg <- readRDS(mgpath)
n <- 200

generate <- function(n, seed) {
    set.seed(seed)
    data_generate(n, tb, mg, reparam = reparam, verbose = FALSE)
}

failures <- 0
seeds <- 1:10
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
  "output/p", p, "q", q, "n", n,
  "reparam", reparam, "sparse", sg, "trials", ntrials, ".rds")

cat("Output will be saved to", path, "\n")

gmmreg_result_mx <- matrix(nrow = ntrials, ncol = 7)
cggr_result_mx <- matrix(nrow = ntrials, ncol = 7)

tictoc::tic()
for (i in seq_along(s_list)) {
  g_result <- gmmreg(
    s_list[[i]]$X, s_list[[i]]$U, 0.75,
    verbose = FALSE, parallel = TRUE)
  gmmreg_result_mx[i, ] <- unlist(performance(g_result, s_list[[i]]))

  c_result <- cggr(
    s_list[[i]]$X, s_list[[i]]$U, 0.75,
    nregmean = 20,
    verbose = FALSE, parallel = TRUE)
  cggr_result_mx[i, ] <- unlist(performance(c_result, s_list[[i]]))

  saveRDS(
    list(
      gmmreg = gmmreg_result_mx,
      cggr = cggr_result_mx,
      run = list(s = s_list[[i]], g_result = g_result, c_result = c_result),
      seeds = seeds),
    path)
  cat("Finished round", i, "\n")

  rm(c_result)
  rm(g_result)
  gc()

}
cat("Finished analysis\n")
tictoc::toc()
cat("Saved results to", path, fill = TRUE)

avg_perf <- rbind(apply(gmmreg_result_mx, 2, mean),
                  apply(gmmreg_result_mx, 2, sd),
                  apply(cggr_result_mx, 2, mean),
                  apply(cggr_result_mx, 2, sd))
colnames(avg_perf) <- c(
  "betaTPR", "betaFPR", "beta_err",
  "meanTPR", "meanFPR", "mean_err", "omega_err")
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
