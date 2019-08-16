context("mpm_to_")

test_that("mpm_to_ functions work correctly", {
  
  xmax <- 20
  
  # mpm_to_lx
  lx <- mpm_to_lx(mat_u, start = 1, xmax = xmax, lx_crit = 0)
  expect_length(lx, xmax + 1)
  expect_true(all(lx >= 0))
  expect_true(all(lx <= 1))
  expect_true(all(lx == cummin(lx))) # monotonic declining
  
  lx_zero <- suppressWarnings(mpm_to_lx(mat_u_zero, start = 1, xmax = xmax))
  expect_equal(lx_zero[1], 1)
  expect_true(all(lx_zero[-1] == 0))
  
  # specify start as integer vs. vector
  lx_a <- mpm_to_lx(mat_u, start = 2)
  lx_b <- mpm_to_lx(mat_u, start = c(0, 1, 0, 0))
  expect_equal(lx_a, lx_b)
  
  # mpm_to_px
  px <- mpm_to_px(mat_u, start = 1, xmax = xmax)
  expect_length(px, length(lx))
  
  # mpm_to_hx
  hx <- mpm_to_hx(mat_u, start = 1, xmax = xmax)
  expect_length(hx, length(lx))
  
  # mpm_to_mx
  mx <- mpm_to_mx(mat_u, mat_f, start = 1, xmax = xmax)
  expect_length(mx, length(lx))
  expect_true(all(mx >= 0))
  
  mx_zero <- suppressWarnings(mpm_to_mx(mat_u_zero, mat_f,
                                        start = 1, xmax = xmax))
  expect_true(all(mx_zero == 0))
})


test_that("mpm_to_ functions warn and fail gracefully", {
  
  # mpm_to_lx
  expect_error(mpm_to_lx(mat_u, start = 10))
  expect_error(mpm_to_lx(mat_u_na))
  
  # mpm_to_px
  expect_error(mpm_to_px(mat_u, start = 10))
  expect_error(mpm_to_px(mat_u_na))
  
  # mpm_to_hx
  expect_error(mpm_to_hx(mat_u, start = 10))
  expect_error(mpm_to_hx(mat_u_na))
  
  # mpm_to_mx
  expect_error(mpm_to_mx(mat_u, mat_f, start = 10))
  expect_error(mpm_to_mx(mat_u_na, mat_f))
  expect_error(mpm_to_mx(mat_u, mat_f_na))
})
