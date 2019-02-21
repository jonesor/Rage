context("R0")

test_that("R0 works correctly", {
  
  x <- R0(mat_u, mat_f)
  x_start <- R0(mat_u, mat_f, method = "startLife")
  x_zero <- suppressWarnings(R0(mat_u_zero, mat_f))
  x_singular <- R0(mat_u_singular, mat_f)

  expect_length(x, 1L)
  expect_true(x > 0)
  expect_length(x_start, 1L)
  expect_true(x_start > 0)
  expect_equal(x_zero, 0)
  expect_equal(x_singular, NA_real_)
})

test_that("R0 warns and fails gracefully", {
  
  expect_warning(R0(mat_u_survissue, mat_f))
  expect_error(R0(mat_u_na, mat_f))
  expect_error(R0(mat_u, mat_f_na))
  expect_error(R0(mat_u, mat_f, method = "eigen"))
})
