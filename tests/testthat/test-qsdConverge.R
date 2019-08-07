context("qsd_converge")

test_that("qsd_converge works correctly", {
  
  x <- qsd_converge(mat_u)
  x_zero <- suppressWarnings(qsd_converge(mat_u_zero, ergodicFix = TRUE))
  
  expect_length(x, 1L)
  expect_true(x > 0)
  expect_equal(x_zero, 1L)
})

test_that("qsd_converge warns and fails gracefully", {
  expect_warning(qsd_converge(mat_u_survissue))
  expect_error(qsd_converge(mat_u_na))
  expect_error(qsd_converge(mat_u, startLife = 10))
})
