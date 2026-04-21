A <- matrix(c(
  0,   0,    2.5, 5.0, 3.0, 1.5, 0.5, 0,
  0.9, 0,    0,   0,   0,   0,   0,   0,
  0,   0.8,  0,   0,   0,   0,   0,   0,
  0,   0,    0.75,0,   0,   0,   0,   0,
  0,   0,    0,   0.7, 0,   0,   0,   0,
  0,   0,    0,   0,   0.6, 0,   0,   0,
  0,   0,    0,   0,   0,   0.5, 0,   0,
  0,   0,    0,   0,   0,   0,   0.4, 0
), nrow = 8, byrow = TRUE)

testthat::expect_error(
  leslie_collapse(cbind(A, A), -4)
)

testthat::expect_true(
  inherits(leslie_collapse(A, 3), "list")
)
