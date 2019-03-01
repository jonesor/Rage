context("mpm_to_")

test_that("mpm_to_ functions work correctly", {
  
  N <- 25
  
  # mpm_to_lx
  lx <- mpm_to_lx(mat_u, startLife = 1, N = N)
  expect_length(lx, N+1)
  expect_true(all(lx >= 0))
  expect_true(all(lx <= 1))
  expect_true(all(lx == cummin(lx))) # monotonic declining
  
  lx_zero <- suppressWarnings(mpm_to_lx(mat_u_zero, startLife = 1, N = N))
  expect_equal(lx_zero[1], 1)
  expect_true(all(lx_zero[-1] == 0))
  
  # mpm_to_px
  px <- mpm_to_px(mat_u, startLife = 1, N = N)
  expect_length(px, N+1)
  
  # mpm_to_hx
  hx <- mpm_to_hx(mat_u, startLife = 1, N = N)
  expect_length(hx, N+1)
  
  # mpm_to_mx
  mx <- mpm_to_mx(mat_u, mat_f, startLife = 1, N = N)
  expect_length(mx, N+1)
  expect_true(all(mx >= 0))
  
  mx_zero <- suppressWarnings(mpm_to_mx(mat_u_zero, mat_f,
                                        startLife = 1, N = N))
  expect_true(all(mx_zero == 0))
})


test_that("mpm_to_ functions warn and fail gracefully", {
  
  # mpm_to_lx
  expect_error(mpm_to_lx(mat_u))
  expect_error(mpm_to_lx(mat_u, startLife = 1))
  expect_error(mpm_to_lx(mat_u, startLife = 10, N = 100))
  expect_error(mpm_to_lx(mat_u_na, startLife = 1, N = 100))
  
  # mpm_to_px
  expect_error(mpm_to_px(mat_u))
  expect_error(mpm_to_px(mat_u, startLife = 1))
  expect_error(mpm_to_px(mat_u, startLife = 10, N = 100))
  expect_error(mpm_to_px(mat_u_na, startLife = 1, N = 100))
  
  # mpm_to_hx
  expect_error(mpm_to_hx(mat_u))
  expect_error(mpm_to_hx(mat_u, startLife = 1))
  expect_error(mpm_to_hx(mat_u, startLife = 10, N = 100))
  expect_error(mpm_to_hx(mat_u_na, startLife = 1, N = 100))
  
  # mpm_to_mx
  expect_error(mpm_to_mx(mat_u, mat_f))
  expect_error(mpm_to_mx(mat_u, mat_f, startLife = 1))
  expect_error(mpm_to_mx(mat_u, mat_f, startLife = 10, N = 100))
  expect_error(mpm_to_mx(mat_u_na, mat_f, startLife = 1, N = 100))
  expect_error(mpm_to_mx(mat_u, mat_f_na, startLife = 1, N = 100))
})
