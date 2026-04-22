data(leslie_mpm1)
A <- leslie_mpm1$matU + leslie_mpm1$matF

testthat::expect_error(
  leslie_collapse(cbind(A, A), -4)
)

testthat::expect_true(
  inherits(leslie_collapse(A, 3), "list")
)

# plus-group survival (A[n,n] > 0) must be preserved after expansion
result <- leslie_collapse(A, 3)
testthat::expect_true(result$Ak[nrow(result$Ak), ncol(result$Ak)] > 0)
