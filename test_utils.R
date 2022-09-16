library(tibble)
library(ggplot2)
library(patchwork)

gamma_viz <- function(gamma_mx, title = "", limits = NULL) {
  gamma_tbl <- as_tibble(cbind(expand.grid(rev(seq_len(d)), seq_len(p)),
                               c(gamma_mx))) |>
               setNames(c("row", "col", "value"))

  ggplot(gamma_tbl, mapping = aes(x = col, y = row, fill = value)) +
    geom_tile(color = "gray30") +
    scale_fill_gradient2(limits = limits) +
    coord_fixed() +
    labs(title = title) +
    theme_minimal()
}

gamma_viz_compare <- function(computed, actual) {
  lim <- max(abs(c(as.numeric(computed), as.numeric(actual))))
  computed_plot <- gamma_viz(computed, "Computed gamma", c(-lim, lim)) +
    guides(fill = "none")
  actual_plot <- gamma_viz(actual, "Actual gamma", c(-lim, lim))
  computed_plot + actual_plot
}

beta_viz <- function(beta_mx, title = "", limits = NULL) {
  beta0 <- as_tibble(cbind(expand.grid(rev(seq_len(d)), seq_len(d)),
                           c(beta_mx))) |>
           setNames(c("row", "col", "value"))

  ggplot(beta0, mapping = aes(x = col, y = row, fill = value)) +
    geom_tile(color = "gray30") +
    scale_fill_gradient2(limits = limits) +
    coord_fixed() +
    labs(title = title) +
    theme_minimal()
}

beta_viz_compare <- function(computed, actual, cov_lbl) {
  lim <- max(abs(c(as.numeric(computed), as.numeric(actual))))
  comp <- beta_viz(computed, title = paste("Computed", cov_lbl),
                   limits = c(-lim, lim))
  act <- beta_viz(actual, title = paste("Actual", cov_lbl),
                  limits = c(-lim, lim))
  comp + act
}

cv_result_plot <- function(cv_result) {
  df <- expand.grid(alpha = seq_along(cv_result$alpha_seq),
                    lambda = seq_along(cv_result$lambda_seq))
  df$error <- as.numeric(cv_result$cv_error)
  opt_lambda_idx <- which(cv_result$lambda_seq == cv_result$opt$lambda)
  opt_alpha_idx <- which(cv_result$alpha == cv_result$opt$alpha)
  ggplot(df) +
    geom_tile(color = "gray80",
              mapping = aes(x = alpha, y = lambda, fill = error)) +
    geom_point(mapping = aes(x = x, y = y),
               data = data.frame(x = opt_alpha_idx, y = opt_lambda_idx),
               color = "tomato") +
    scale_x_continuous(breaks = seq_along(cv_result$alpha_seq),
                       labels = cv_result$alpha_seq) +
    scale_y_reverse(breaks = seq_along(cv_result$lambda_seq),
                       labels = round(cv_result$lambda_seq, 3)) +
    scale_fill_gradient(low = "white", high = "black") +
    coord_fixed() +
    theme_classic()
}

#' @return true rate of zero coefficients excluding diagonals
true_neg_rate <- function(est_beta_mx, true_beta_mx, true_cov_nz_idx) {
  stopifnot(is.list(est_beta_mx),
            is.list(true_beta_mx),
            length(true_cov_nz_idx) + 1 == length(true_beta_mx))

  # Counts the number of entries that are simultaneously zero in both a and b,
  # excluding diagonal entries. If b is NULL, counts the number of off-diag
  # zeros in a.
  off_diag_zero <- function(a, b = NULL) {
    if (is.null(b)) {
      b <- matrix(0, dim(a)[1], dim(a)[2])
    }
    stopifnot(dim(a) == dim(b))
    diag(a) <- NA
    diag(b) <- NA
    sum(a == 0 & b == 0, na.rm = T)
  }

  # Index of non-zero groups, including the population group (hence + 1)
  true_nz_groups <- c(1, 1 + true_cov_nz_idx)
  matching_nz_groups <- map2(est_beta_mx[true_nz_groups], true_beta_mx,
                         ~off_diag_zero(.x, .y)) |>
                    unlist() |> sum()
  matching_z_groups <- map(est_beta_mx[-true_nz_groups], ~off_diag_zero(.x)) |>
                unlist() |> sum()
  true_test_zero <- matching_nz_groups + matching_z_groups

  d <- nrow(true_beta_mx[[1]])
  true_zero <- sum(unlist(map(true_beta_mx, ~off_diag_zero(.x)))) +
    (length(est_beta_mx) - 1 - length(true_cov_nz_idx)) * (d^2 - d)

  result <- true_test_zero / true_zero
  return(round(result, 4))
}

#' @return true rate of non-zero coefficients excluding diagonals
true_pos_rate <- function(est_beta_mx, true_beta_mx, true_cov_nz_idx) {
  stopifnot(is.list(est_beta_mx),
            is.list(true_beta_mx),
            length(true_cov_nz_idx) + 1 == length(true_beta_mx))

  # Counts the number of entries that are simultaneously nonzero in both a and b,
  # excluding diagonal entries. If b is NULL, counts the number of off-diag
  # nonzeros in a.
  off_diag_nonzero <- function(a, b = NULL) {
    if (is.null(b)) {
      b <- matrix(1, dim(a)[1], dim(a)[2])
    }
    stopifnot(dim(a) == dim(b))
    diag(a) <- NA
    diag(b) <- NA
    sum(a != 0 & b != 0, na.rm = T)
  }

  # Index of non-zero groups, including the population group (hence + 1)
  true_nz_groups <- c(1, 1 + true_cov_nz_idx)
  matching_nz_groups <- map2(est_beta_mx[true_nz_groups], true_beta_mx,
                             ~off_diag_nonzero(.x, .y)) |>
    unlist() |> sum()
  true_test_nonzero <- matching_nz_groups

  true_nonzero <- sum(unlist(map(true_beta_mx, ~off_diag_nonzero(.x))))

  result <- true_test_nonzero / true_nonzero
  return(round(result, 4))
}
