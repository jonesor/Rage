data(leslie_mpm1)
A <- leslie_mpm1$matU + leslie_mpm1$matF

testthat::expect_error(
Leslie_collapse(cbind(A,A),-4)
)

testthat::expect_true(
inherits(Leslie_collapse(A,3),"list")
)
