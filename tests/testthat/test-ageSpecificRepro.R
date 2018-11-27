context("ageSpecificRepro")

test_that("ageSpecificRepro works correctly", {
  
  N <- 100
  
  x <- ageSpecificRepro(mat_u, mat_f, startLife = 1, N = N)
  x_zero <- suppressWarnings(
    ageSpecificRepro(mat_u_zero, mat_f, startLife = 1, N = N)
  )
  
  expect_length(x, N+1)
  expect_true(all(x >= 0))
  
  expect_true(all(x_zero == 0))
})

test_that("ageSpecificRepro warns and fails gracefully", {
  
  expect_error(ageSpecificRepro(mat_u, mat_f))
  expect_error(ageSpecificRepro(mat_u, mat_f, startLife = 1))
  expect_error(ageSpecificRepro(mat_u, mat_f, startLife = 10, N = 100))
  expect_error(ageSpecificRepro(mat_u_na, mat_f, startLife = 1, N = 100))
  expect_error(ageSpecificRepro(mat_u, mat_f_na, startLife = 1, N = 100))
})
