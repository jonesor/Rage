test_that("mpm_to_ functions work correctly", {
  xmax <- 20

  # mpm_to_lx
  lx <- mpm_to_lx(mat_u, start = 1, xmax = xmax, lx_crit = 0)
  expect_length(lx, xmax + 1)
  expect_true(all(lx >= 0))
  expect_true(all(lx <= 1))
  expect_true(all(lx == cummin(lx))) # monotonic declining

  lx_zero <- suppressWarnings(mpm_to_lx(mat_u_zero, start = 1, xmax = xmax))
  expect_identical(lx_zero[1], 1)
  expect_true(all(lx_zero[-1] == 0))

  # specify start as integer vs. vector
  lx_a <- mpm_to_lx(mat_u, start = 2)
  lx_b <- mpm_to_lx(mat_u, start = c(0, 1, 0, 0))
  expect_identical(lx_a, lx_b)

  # mpm_to_px
  px <- mpm_to_px(mat_u, start = 1, xmax = xmax, lx_crit = 0)
  expect_length(px, length(lx))

  # mpm_to_hx
  hx <- mpm_to_hx(mat_u, start = 1, xmax = xmax, lx_crit = 0)
  expect_length(hx, length(lx))

  # mpm_to_mx
  mx <- mpm_to_mx(mat_u, mat_f, start = 1, xmax = xmax, lx_crit = 0)
  expect_length(mx, length(lx))
  expect_true(all(mx >= 0))

  mx_zero <- suppressWarnings(mpm_to_mx(mat_u_zero, mat_f,
    start = 1, xmax = xmax
  ))
  expect_true(all(mx_zero == 0))

  # using named life stages
  lx_named <- mpm_to_lx(mat_u_named, start = "sm", xmax = xmax, lx_crit = 0)
  px_named <- mpm_to_px(mat_u_named, start = 1, xmax = xmax, lx_crit = 0)
  hx_named <- mpm_to_hx(mat_u_named, start = 1, xmax = xmax, lx_crit = 0)
  mx_named <- mpm_to_mx(mat_u_named, mat_f_named,
    start = 1, xmax = xmax,
    lx_crit = 0
  )
  expect_identical(lx, lx_named)
  expect_identical(px, px_named)
  expect_identical(hx, hx_named)
  expect_identical(mx, mx_named)
})


test_that("mpm_to_ functions warn and fail gracefully", {
  # mpm_to_lx
  expect_error(mpm_to_lx(mat_u, start = 10))
  expect_error(mpm_to_lx(mat_u_na))
  expect_error(mpm_to_lx(mat_u, start = "stage name"))

  # mpm_to_px
  expect_error(mpm_to_px(mat_u, start = 10))
  expect_error(mpm_to_px(mat_u_na))
  expect_error(mpm_to_px(mat_u, start = "stage name"))

  # mpm_to_hx
  expect_error(mpm_to_hx(mat_u, start = 10))
  expect_error(mpm_to_hx(mat_u_na))
  expect_error(mpm_to_hx(mat_u, start = "stage name"))

  # mpm_to_mx
  expect_error(mpm_to_mx(mat_u, mat_f, start = 10))
  expect_error(mpm_to_mx(mat_u_na, mat_f))
  expect_error(mpm_to_mx(mat_u, mat_f_na))
  expect_error(mpm_to_mx(mat_u, mat_f, start = "stage name"))
})
