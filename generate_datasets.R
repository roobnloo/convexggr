source("data_generation.R")
source("cggr.R")
source("gmmreg.R")
source("test_utils.R")

# generate_parameters(25, 100, 10, FALSE)

p <- 25
q <- 100
tb <- readRDS("data/tb_p25q100.rds")
mg <- readRDS("data/mg_p25q100sparseFALSE.rds")
n <- 400

generate <- function(n, seed) {
    set.seed(seed)
    s <- data_generate(n, tb, mg, reparam = TRUE, verbose = FALSE)
}

failures <- 0
trials <- 1
s_list <- vector(mode = "list", length = trials)
for (i in seq_len(trials)) {
  s <- generate(n, i)
  if (is.null(s)) {
    failures <- failures + 1
    next
  }
  s_list[[i]] <- s
}

cat("Done generating datasets with", failures, "failures\n")

gmmreg_result_mx <- matrix(nrow = trials, ncol = 5)
cggr_result_mx <- matrix(nrow = trials, ncol = 5)
tictoc::tic()
for (i in seq_along(s_list)) {
  # g_result <- gmmreg(
  #   s_list[[i]]$X, s_list[[i]]$U, 0.75,
  #   verbose = FALSE, parallel = TRUE)
  # gmmreg_result_mx[i, ] <- unlist(performance(g_result, s_list[[i]]))

  c_result <- cggr(
    s_list[[i]]$X, s_list[[i]]$U, 0.75, regmean = 3.00,
    verbose = TRUE, parallel = TRUE)
  cggr_result_mx[i, ] <- unlist(performance(c_result, s_list[[i]]))
  cat("Finished round", i, "\n")
}
cat("Finished analysis\n")
tictoc::toc()

results <- list(
    gmmreg = gmmreg_result_mx,
    cggr = cggr_result_mx)
# path <- paste0("output/p", p, "q", q, "n", n, ".rds")
# saveRDS(results, path)
# cat("Saved results to", path, "\n")

avg_perf <- rbind(apply(gmmreg_result_mx, 2, mean),
                  apply(gmmreg_result_mx, 2, sd),
                  apply(cggr_result_mx, 2, mean),
                  apply(cggr_result_mx, 2, sd))
colnames(avg_perf) <- c("TPR", "FPR", "beta_err", "mean_err", "omega_err")
rownames(avg_perf) <- c("gmmreg", "(sd)", "cggr", "(sd)")

# for (i in seq_len(6)) {
#   cov_lbl <- as.character(i - 1)
#   beta_list <- list(GMMReg = g_result$bhat[, , i],
#                     CGGR = c_result$bhat[, , i],
#                     Actual = tb[, , i])
#   print(beta_viz_list(beta_list, cov_lbl, tileborder = TRUE))
# }

# print(gamma_viz_list(list(g = g_result$ghat, c = c_result$ghat, mg)))
