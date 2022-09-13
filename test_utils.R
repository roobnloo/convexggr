library(tibble)
library(ggplot2)

gamma_viz <- function(gamma_mx) {
  gamma_tbl <- as_tibble(cbind(expand.grid(rev(seq_len(d)), seq_len(p)),
                               c(gamma_mx))) |>
               setNames(c("row", "col", "value"))

  ggplot(gamma_tbl, mapping = aes(x = col, y = row, fill = value)) +
    geom_tile(color = "gray30") +
    scale_fill_gradient2() +
    coord_fixed() +
    theme_minimal()
}

beta_viz <- function(beta_mx, title = "") {
  beta0 <- as_tibble(cbind(expand.grid(rev(seq_len(d)), seq_len(d)),
                           c(beta_mx))) |>
           setNames(c("row", "col", "value"))

  ggplot(beta0, mapping = aes(x = col, y = row, fill = value)) +
    geom_tile(color = "gray30") +
    scale_fill_gradient2() +
    coord_fixed() +
    labs(title = title) +
    theme_minimal()
}

#' @return  the True Positive Rate of estimated zeros compared to the actual
#' support set. Diagonal entries are ignored.
true_pos_rate <- function(est_beta_mx, true_beta_mx, true_cov_nz_idx) {
  stopifnot(is.list(est_beta_mx),
            is.list(true_beta_mx),
            length(true_cov_nz_idx) + 1 == length(true_beta_mx))

  off_diag_zero <- function(a, b = NULL) {
    if (is.null(b)) {
      b <- matrix(0, dim(a)[1], dim(a)[2])
    }
    stopifnot(dim(a) == dim(b))
    diag(a) <- NA
    diag(b) <- NA
    sum(a == 0 & b == 0, na.rm = T)
  }

  matching_idx <- c(1, 1 + true_cov_nz_idx)
  true_test_zero <- map2(est_beta_mx[matching_idx], true_beta_mx,
                         ~off_diag_zero(.x, .y)) |>
                    unlist() |> sum()
  additional <- map(est_beta_mx[-matching_idx], ~off_diag_zero(.x)) |>
                unlist() |> sum()
  true_test_zero <- true_test_zero + additional

  d <- nrow(true_beta_mx[[1]])
  p <- length(est_beta_mx) - 1
  true_zero <- sum(unlist(map(true_beta_mx, ~off_diag_zero(.x)))) +
    (p - length(true_cov_nz_idx)) * (d^2 - d)

  result <- true_test_zero / true_zero
  return(result)
}
