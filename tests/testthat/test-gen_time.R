
test_that("gen_time works correctly", {
  
  x <- gen_time(mat_u, mat_f)
  x_start <- gen_time(mat_u, mat_f)
  x_singular <- gen_time(mat_u_singular, mat_f)

  expect_length(x, 1L)
  expect_true(x > 0)
  expect_length(x_start, 1L)
  expect_true(x_start > 0)
  expect_equal(x_singular, NA_real_)
  
  x_po <- gen_time(mat_u, mat_f, method = "parent_offspring")
  x_singular_po <- gen_time(mat_u_singular, mat_f, method = "parent_offspring")
  expect_length(x_po, 1L)
  expect_true(x_po > 0)
  expect_equal(x_singular_po, NA_real_)
  
  x_co <- gen_time(mat_u, mat_f, method = "cohort")
  x_singular_co <- gen_time(mat_u_singular, mat_f, method = "cohort")
  expect_length(x_co, 1L)
  expect_true(x_co > 0)
  expect_equal(x_singular_co, NA_real_)
})

test_that("gen_time warns and fails gracefully", {
  
  expect_warning(gen_time(mat_u_survissue, mat_f))
  expect_error(gen_time(mat_u_na, mat_f))
  expect_error(gen_time(mat_u, mat_f_na))
  
  expect_warning(gen_time(mat_u_survissue, mat_f), method = "parent_offspring")
  expect_error(gen_time(mat_u_na, mat_f), method = "parent_offspring")
  expect_error(gen_time(mat_u, mat_f_na), method = "parent_offspring")
  
  expect_warning(gen_time(mat_u_survissue, mat_f), method = "cohort")
  expect_error(gen_time(mat_u_na, mat_f), method = "cohort")
  expect_error(gen_time(mat_u, mat_f_na), method = "cohort")
})
