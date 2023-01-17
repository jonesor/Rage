
test_that("net_repro_rate works correctly", {
  x <- net_repro_rate(mat_u, mat_f)
  x_start <- net_repro_rate(mat_u, mat_f, method = "start")
  x_zero <- suppressWarnings(net_repro_rate(mat_u_zero, mat_f))
  x_singular <- net_repro_rate(mat_u_singular, mat_f)
  x_named <- net_repro_rate(mat_u_named, mat_f_named, start = "sm")

  expect_length(x, 1L)
  expect_true(x > 0)
  expect_length(x_start, 1L)
  expect_true(x_start > 0)
  expect_equal(x_zero, 0)
  expect_equal(x_singular, NA_real_)
  expect_equal(x, x_named)
})

test_that("net_repro_rate warns and fails gracefully", {
  expect_warning(net_repro_rate(mat_u_survissue, mat_f))
  expect_error(net_repro_rate(mat_u_na, mat_f))
  expect_error(net_repro_rate(mat_u, mat_f_na))
  expect_error(net_repro_rate(mat_u, mat_f, method = "eigen"))
  expect_error(net_repro_rate(mat_u_named, mat_f, start = "sm"))
  expect_error(net_repro_rate(mat_u, mat_f_named, start = "sm"))
  expect_error(net_repro_rate(mat_u_named, mat_f_named, start = "invalid stage"))
})
