context("entropy_d")

test_that("entropy_d works correctly", {

  lx <- 0.8^(0:8)
  mx <- c(0.00, 0.00, 1.10, 1.50, 1.60, 1.70, 1.50, 1.20, 0.70)
  
  x <- entropy_d(lx, mx)

  expect_length(x, 1L)
  expect_true(x > 0)
})

test_that("entropy_d warns and fails gracefully", {
  
  lx1 <- c(1.1, 0.6, 0.5, 0.4)  # lx > 1
  mx1 <- c(4.4, 3.3, 2.2, 1.1)
  expect_error(entropy_d(lx1, mx1))
  
  lx2 <- c(1.0, 0.6, 0.61, 0.5) # lx not monotonically declining
  mx2 <- c(4.4, 3.3, 2.2, 1.1)
  expect_error(entropy_d(lx2, mx2))
  
  lx3 <- c(1.0, 0.6, 0.6, 0.5) # negative mx
  mx3 <- c(4.4, 3.3, 2.2, -0.1)
  expect_error(entropy_d(lx3, mx3))
})
