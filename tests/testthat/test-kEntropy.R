context("kEntropy")

test_that("kEntropy works correctly", {
  
  x <- kEntropy(mat_u)
  x_trap <- kEntropy(mat_u, trapeze = TRUE)

  expect_length(x, 1L)
  expect_length(x_trap, 1L)
  expect_true(x >= 0)
  expect_true(x_trap >= 0)
})

test_that("kEntropy warns and fails gracefully", {
  
  expect_warning(kEntropy(mat_u_survissue))
  expect_error(kEntropy(mat_u_na))
  expect_error(kEntropy(mat_u, startLife = 10))
})
