context("qsdConverge")

test_that("qsdConverge works correctly", {
  
  x <- qsdConverge(mat_u)
  x_zero <- suppressWarnings(qsdConverge(mat_u_zero))
  
  expect_length(x, 1L)
  expect_true(x > 0)
  expect_equal(x_zero, NA_integer_)
})

test_that("qsdConverge warns and fails gracefully", {
  
  expect_warning(qsdConverge(mat_u_survissue))
  expect_error(qsdConverge(mat_u_na))
  expect_error(qsdConverge(mat_u, startLife = 10))
})
