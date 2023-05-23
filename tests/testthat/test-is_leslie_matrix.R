A <- matrix(c(
  0.1, 1.2, 1.1,
  0.1, 0.0, 0.0,
  0.0, 0.2, 0.3
), nrow = 3, byrow = TRUE)

expect_true(is_leslie_matrix(A)) # true

A <- matrix(c(
  0.1, 1.2, 1.1,
  0.1, 0.2, 0.1,
  0.2, 0.3, 0.3
), nrow = 3, byrow = TRUE)

expect_false(is_leslie_matrix(A))
data(leslie_mpm1)

A <- leslie_mpm1$matU + leslie_mpm1$matF

expect_true(is_leslie_matrix(A)) # false
