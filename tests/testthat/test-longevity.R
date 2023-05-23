test_that("longevity works correctly", {
  x <- longevity(mat_u)
  x_zero <- suppressWarnings(longevity(mat_u_zero))
  x_noconverg <- suppressWarnings(longevity(mat_u_singular, start = 3))
  x_named <- longevity(mat_u_named, start = "sm")

  expect_length(x, 1L)
  expect_gt(x, 0)
  expect_identical(x_zero, 1L)
  expect_identical(as.integer(x_noconverg), as.integer(NA_real_))
  expect_identical(x, x_named)
})

test_that("longevity warns and fails gracefully", {
  expect_warning(longevity(mat_u_zero))
  expect_warning(longevity(mat_u_survissue))
  expect_warning(longevity(mat_u_singular, start = 3))
  expect_error(longevity(mat_u_na))
  expect_error(longevity(mat_u, start = 10))
  expect_error(longevity(mat_u, start = "stage name"))
  expect_error(longevity(mat_u_named, start = "invalid stage"))
  expect_error(longevity(mat_u, lx_crit = 10))
})
