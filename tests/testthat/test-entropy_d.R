context("entropy_d")

test_that("entropy_d works correctly", {
  
  x <- entropy_d(mat_u, mat_f)

  expect_length(x, 1L)
  expect_true(x > 0)
})

test_that("entropy_d warns and fails gracefully", {
  
  expect_warning(entropy_d(mat_u_survissue, mat_f))
  expect_error(entropy_d(mat_u_na, mat_f))
  expect_error(entropy_d(mat_u, mat_f_na))
  expect_error(entropy_d(mat_u, mat_f, startLife = 10))
})
