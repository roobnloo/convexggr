library(assertthat)

x0 <- matrix(1, 2, 2)
y0 <- matrix(1, 2, 2)

x1 <- matrix(0, 2, 2)
y1 <- matrix(c(0, 0, 1, 1), 2, 2)

assert_that(true_pos_rate(list(x0, x1), list(y0, y1), 1) == 1)

assert_that(true_pos_rate(list(x0, x1), list(diag(2), y1), 1) == 1/3)

assert_that(true_pos_rate(list(x0, x0), list(diag(2), y1), 1) == 0)

x2 <- matrix(0:1, 2, 2)
x3 <- matrix(1, 2, 2)

assert_that(true_pos_rate(list(x0, x1, x2, x3),
                          list(x0, x1, x2, x3), 1:3) == 1)
assert_that(true_pos_rate(list(x0, x1, x2, x3),
                          list(x0, x1, x2), 1:2) == 3/5)
assert_that(true_pos_rate(list(x0, x1, x2, x3),
                          list(x0, x1), 1) == 3/6)


map(seq_along(result), ~ apply(abs(result[[.x]]), 1, sum)) |>
  reduce(`+`)
