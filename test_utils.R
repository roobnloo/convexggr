library(tibble)
library(ggplot2)
library(patchwork)

performance <- function(result, s) {
  stopifnot(all(dim(result$bhat) == dim(s$tB)))
  p <- dim(s$tB)[1]
  n <- length(s$pd_check)
  stats <- list()
  stats$tpr <- sum(result$bhat != 0 & s$tB != 0) / sum(s$tB != 0)
  stats$fpr <- sum(result$bhat != 0 & s$tB == 0) / sum(s$tB == 0)

  beta_err <- 0
  for (i in 1:p) {
    beta_err <- beta_err + sqrt(sum((s$tB[i, , ] - result$bhat[i, , ])^2))
  }
  stats$beta_err <- beta_err

  mean_err <- 0
  for (i in 1:n) {
    mean_err <- mean_err + sum((s$mumx[i, ] -  result$mean(i))^2) / n
  }
  stats$mean_err <- mean_err / n


  omega_err <- 0
  iu <- cbind(1, s$U)
  for (i in 1:n) {
    omega <- -apply(s$tB, c(1, 2), \(b) b %*% iu[i, ])
    omhat <- result$precision(i)
    diag(omhat) <- 0
    omega_err <- omega_err + sum((omega - omhat)^2) / n
  }
  stats$omega_err <- omega_err

  return(stats)
}

gamma_viz <- function(gamma_mx, title = "", limits = NULL) {
  d <- nrow(gamma_mx)
  p <- ncol(gamma_mx)
  gamma_tbl <- as_tibble(cbind(expand.grid(rev(seq_len(d)), seq_len(p)),
                               c(gamma_mx))) |>
               setNames(c("row", "col", "value"))

  ggplot(gamma_tbl, mapping = aes(x = col, y = row, fill = value)) +
    # geom_tile(color = "gray30") +
    geom_tile() +
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

gamma_viz_list <- function(gamma_list) {
  lim <- max(abs(unlist(gamma_list)))
  plots <- vector(mode = "list", length = length(gamma_list))
  for (i in seq_along(gamma_list)) {
    plots[[i]] <- gamma_viz(gamma_list[[i]],
                            paste(names(gamma_list)[i], "gamma"), c(-lim, lim)) +
                  labs(x = NULL, y = NULL) +
                  guides(x = "none", y = "none")
    if (i != length(gamma_list)) {
      plots[[i]] <- plots[[i]] + guides(fill = "none")
    }
  }
  Reduce(`+`, plots) + plot_layout(ncol = 1)
}

beta_viz <- function(beta_mx, title = "", limits = NULL, guides = T,
                     tileborder = T, fill_legend = T) {
  d <- nrow(beta_mx)
  beta0 <- as_tibble(cbind(expand.grid(rev(seq_len(d)), seq_len(d)),
                           c(beta_mx))) |>
           setNames(c("row", "col", "value"))

  tilecolor <- ifelse(tileborder, "gray30", "white")
  p <- ggplot(beta0, mapping = aes(x = col, y = row, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(limits = limits) +
        coord_fixed() +
        labs(title = title) +
        theme_minimal()
  if (!fill_legend) {
    p <- p + guides(fill = "none")
  }
  else if (!guides) {
    p <- p + guides(x = "none", y = "none", fill = "none") +
          labs(x = NULL, y = NULL)
  }

  return(p)
}

beta_viz_compare <- function(computed, actual, cov_lbl, guides = T,
                             tileborder = T) {
  lim <- max(abs(c(as.numeric(computed), as.numeric(actual))))
  comp <- beta_viz(computed, title = paste("Computed", cov_lbl),
                   limits = c(-lim, lim),
                   guides = guides,
                   tileborder = tileborder,
                   fill_legend = F)
  act <- beta_viz(actual, title = paste("Actual", cov_lbl),
                  limits = c(-lim, lim),
                  guides = guides,
                  tileborder = tileborder)
  comp + act
}

beta_viz_list <- function(beta_list, cov_lbl, guides = T, tileborder = T) {
  lim <- max(abs(unlist(beta_list)))

  plots <- vector(mode = "list", length = length(beta_list))
  for (i in seq_along(beta_list)) {
    plots[[i]] <-  beta_viz(beta_list[[i]],
                            title = paste(names(beta_list)[i], cov_lbl),
                            limits = c(-lim, lim),
                            tileborder = tileborder) +
                   guides(x = "none", y = "none") +
                   labs(x = NULL, y = NULL)
    if (i != length(beta_list)) {
      plots[[i]] <- plots[[i]] + guides(fill = "none")
    }
  }

  Reduce(`+`, plots)
}

cv_result_plot <- function(cv_result, guides = T) {
  df <- expand.grid(alpha = seq_along(cv_result$alpha_seq),
                    lambda = seq_along(cv_result$lambda_seq))
  df$error <- as.numeric(cv_result$cv_error)
  opt_lambda_idx <- which(cv_result$lambda_seq == cv_result$opt$lambda)
  opt_alpha_idx <- which(cv_result$alpha == cv_result$opt$alpha)
  p <- ggplot(df) +
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
  if (!guides) {
    p <- p + guides(x = "none", y = "none", fill = "none")
  }
  return(p)
}

#' @return true rate of zero coefficients excluding diagonals
true_neg_rate <- function(est_beta_mx, true_beta_mx, qe) {
  stopifnot(is.list(est_beta_mx),
            is.list(true_beta_mx),
            qe + 1 == length(true_beta_mx))

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
  true_nz_groups <- seq_len(qe+1)
  matching_nz_groups <- map2(est_beta_mx[true_nz_groups], true_beta_mx,
                         ~off_diag_zero(.x, .y)) |>
                    unlist() |> sum()
  matching_z_groups <- map(est_beta_mx[-true_nz_groups], ~off_diag_zero(.x)) |>
                unlist() |> sum()
  true_test_zero <- matching_nz_groups + matching_z_groups

  d <- nrow(true_beta_mx[[1]])
  true_zero <- sum(unlist(map(true_beta_mx, ~off_diag_zero(.x)))) +
    (length(est_beta_mx) - 1 - qe) * (d^2 - d)

  result <- true_test_zero / true_zero
  return(round(result, 4))
}

#' @return true rate of non-zero coefficients excluding diagonals
true_pos_rate <- function(est_beta_mx, true_beta_mx, qe) {
  stopifnot(is.list(est_beta_mx),
            is.list(true_beta_mx),
            qe + 1 == length(true_beta_mx))

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
  true_nz_groups <- seq_len(qe + 1)
  matching_nz_groups <- map2(est_beta_mx[true_nz_groups], true_beta_mx,
                             ~off_diag_nonzero(.x, .y)) |>
    unlist() |> sum()
  true_test_nonzero <- matching_nz_groups

  true_nonzero <- sum(unlist(map(true_beta_mx, ~off_diag_nonzero(.x))))

  result <- true_test_nonzero / true_nonzero
  return(round(result, 4))
}
