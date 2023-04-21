data(leslie_mpm1)
A <- leslie_mpm1$matU + leslie_mpm1$matF

testthat::expect_error(
leslie_collapse(cbind(A,A),-4)
)

testthat::expect_true(
inherits(leslie_collapse(A,3),"list")
)