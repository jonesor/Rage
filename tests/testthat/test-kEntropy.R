context("entropy_k")

test_that("entropy_k works correctly", {
  
  x <- entropy_k(mat_u)
  x_trap <- entropy_k(mat_u, trapeze = TRUE)

  expect_length(x, 1L)
  expect_length(x_trap, 1L)
  expect_true(x >= 0)
  expect_true(x_trap >= 0)
})

test_that("entropy_k warns and fails gracefully", {
  
  expect_warning(entropy_k(mat_u_survissue))
  expect_error(entropy_k(mat_u_na))
  expect_error(entropy_k(mat_u, startLife = 10))
})
