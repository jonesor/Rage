
test_that("gen_time works correctly", {
  
  x <- gen_time(mat_u, mat_f)
  x_start <- gen_time(mat_u, mat_f)
  x_singular <- gen_time(mat_u_singular, mat_f)

  expect_length(x, 1L)
  expect_true(x > 0)
  expect_length(x_start, 1L)
  expect_true(x_start > 0)
  expect_equal(x_singular, NA_real_)
})

test_that("gen_time warns and fails gracefully", {
  
  expect_warning(gen_time(mat_u_survissue, mat_f))
  expect_error(gen_time(mat_u_na, mat_f))
  expect_error(gen_time(mat_u, mat_f_na))
})
