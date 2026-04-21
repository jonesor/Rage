A <- matrix(c(
  0.1, 1.2, 1.1,
  0.1, 0.0, 0.0,
  0.0, 0.2, 0.0
), nrow = 3, byrow = TRUE)

expect_true(is_leslie_matrix(A))

# plus-group matrix (A[n,n] > 0) is not a Leslie matrix
A_plus <- matrix(c(
  0.1, 1.2, 1.1,
  0.1, 0.0, 0.0,
  0.0, 0.2, 0.3
), nrow = 3, byrow = TRUE)

expect_false(is_leslie_matrix(A_plus))

A <- matrix(c(
  0.1, 1.2, 1.1,
  0.1, 0.2, 0.1,
  0.2, 0.3, 0.3
), nrow = 3, byrow = TRUE)

expect_false(is_leslie_matrix(A))
data(leslie_mpm1)

A <- leslie_mpm1$matU + leslie_mpm1$matF

# leslie_mpm1 has A[n,n] > 0 (plus-group variant), so is not a Leslie matrix
expect_false(is_leslie_matrix(A))
