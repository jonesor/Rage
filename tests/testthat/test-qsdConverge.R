context("qsdConverge")

test_that("qsdConverge works correctly", {
  
  x <- qsdConverge(mat_u)
  x_zero <- suppressWarnings(qsdConverge(mat_u_zero, ergodicFix = TRUE))
  
  expect_length(x, 1L)
  expect_true(x > 0)
  expect_equal(x_zero, 1L)
})

test_that("qsdConverge warns and fails gracefully", {
  expect_warning(qsdConverge(mat_u_survissue))
  expect_error(qsdConverge(mat_u_na))
  expect_error(qsdConverge(mat_u, startLife = 10))
  expect_error(qsdConverge(mat_u_zero))
})
