source("cv_cggr.R")
source("cggr.R")
library(tictoc)

run <- function(data, alpha_seq, lambda_len, lambda_g = 1,
                out_path = "output/run.out", verbose = T) {
  d <- ncol(data$responses)
  p <- ncol(data$covariates)
  num <- nrow(data$responses)

  tuned_params <- matrix(NA, nrow = 2, ncol = d)
  rownames(tuned_params) <- c("alpha", "lambda")
  cv_results <- vector(mode = "list", length = d)

  if (verbose) {
    print(paste("Cross-validating over a", length(alpha_seq), "x",
                lambda_len, "grid..."))
    tic("Cross-validation")
  }
  for (i in seq_len(d)) {
    cv_result <- cv_cggr(data$responses[, i],
                         data$responses[, -i],
                         data$covariates,
                         lambda_g, nfold = 5, alpha_seq = alpha_seq,
                         lambda_len = lambda_len, lambda_min_prop = 0.1)
    cv_results[[i]] <- cv_result
    tuned_params[, i] <- unlist(cv_result$opt)
  }

  if (verbose) {
    toc()
    print("Running cggr on tuned parameters...")
    tic("cggr")
  }

  cggr_result <- convex_ggr(data$responses, data$covariates, lambda_g,
                            tuned_params[2, ], tuned_params[1, ])

  if(verbose) {
    toc()
  }

  run_result <- list(cv_result = cv_results,
                     tuned = tuned_params,
                     cggr_result = cggr_result)

  print(paste("Saving result to", out_path))
  saveRDS(run_result, "output/run.out")

  return(run_result)
}
