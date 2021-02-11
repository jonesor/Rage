context("entropy_k")

test_that("entropy_k works correctly", {
  
  lx <- 0.8^(0:20)
  
  x <- entropy_k(lx)
  x_trap <- entropy_k(lx, trapeze = TRUE)

  expect_length(x, 1L)
  expect_length(x_trap, 1L)
  expect_true(x >= 0)
  expect_true(x_trap >= 0)
})

test_that("entropy_k warns and fails gracefully", {
  
  lx1 <- c(1.1, 0.6, 0.5, 0.4)  # lx > 1
  expect_error(entropy_k(lx1))
  
  lx2 <- c(1.0, 0.6, 0.61, 0.5) # lx not monotonically declining
  expect_error(entropy_k(lx2))
})
