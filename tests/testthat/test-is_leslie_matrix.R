# Standard Leslie matrix
A <- matrix(c(
  0.1, 1.2, 1.1,
  0.1, 0.0, 0.0,
  0.0, 0.2, 0.0
), nrow = 3, byrow = TRUE)

expect_true(is_leslie_matrix(A))

# Plus-group Leslie matrix: A[n,n] > 0 is permitted
A_plus <- matrix(c(
  0.1, 1.2, 1.1,
  0.1, 0.0, 0.0,
  0.0, 0.2, 0.3
), nrow = 3, byrow = TRUE)

expect_true(is_leslie_matrix(A_plus))

# Not a Leslie matrix: non-zero off-diagonal elements
A <- matrix(c(
  0.1, 1.2, 1.1,
  0.1, 0.2, 0.1,
  0.2, 0.3, 0.3
), nrow = 3, byrow = TRUE)

expect_false(is_leslie_matrix(A))

# leslie_mpm1 is a plus-group Leslie matrix (A[n,n] > 0)
data(leslie_mpm1)

A <- leslie_mpm1$matU + leslie_mpm1$matF

expect_true(is_leslie_matrix(A))
