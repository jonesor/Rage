context("ageSpecificSurv")

test_that("ageSpecificSurv works correctly", {
  
  N <- 100
  
  x <- ageSpecificSurv(mat_u, startLife = 1, N = N)
  expect_length(x, N+1)
  expect_true(all(x >= 0))
  expect_true(all(x <= 1))
  expect_true(all(x == cummin(x))) # monotonic declining
  
  x_zero <- suppressWarnings(ageSpecificSurv(mat_u_zero, startLife = 1, N = N))
  expect_equal(x_zero[1], 1)
  expect_true(all(x_zero[-1] == 0))
})

test_that("ageSpecificSurv warns and fails gracefully", {
  
  expect_error(ageSpecificSurv(mat_u))
  expect_error(ageSpecificSurv(mat_u, startLife = 1))
  expect_error(ageSpecificSurv(mat_u, startLife = 10, N = 100))
  expect_error(ageSpecificSurv(mat_u_na, startLife = 1, N = 100))
})
