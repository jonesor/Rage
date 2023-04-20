test_that("mpm_to_table works correctly", {
  xmax <- 20

  x_u <- mpm_to_table(mat_u, xmax = xmax, lx_crit = 0)
  x_uf <- mpm_to_table(mat_u, mat_f, xmax = xmax, lx_crit = 0)
  x_uc <- mpm_to_table(mat_u, matC = mat_c, xmax = xmax, lx_crit = 0)
  x_ufc <- mpm_to_table(mat_u, mat_f, mat_c, xmax = xmax, lx_crit = 0)
  x_u_named <- mpm_to_table(mat_u_named, start = "sm", xmax = xmax, lx_crit = 0)
  x_uf_named <- mpm_to_table(mat_u_named, mat_f_named,
    start = "sm", xmax = xmax,
    lx_crit = 0
  )

  expect_s3_class(x_u, "data.frame")
  expect_equal(nrow(x_u), xmax + 1)
  expect_equal(ncol(x_u), 7)
  expect_equal(ncol(x_uf), 9)
  expect_equal(ncol(x_uc), 9)
  expect_equal(ncol(x_ufc), 13)
  expect_equal(x_u, x_u_named)
  expect_equal(x_uf, x_uf_named)
})

test_that("mpm_to_table warns and fails gracefully", {
  xmax <- 20

  expect_error(mpm_to_table(mat_u_na, mat_f, xmax = xmax))
  expect_error(mpm_to_table(mat_u, mat_f_na, xmax = xmax))
  expect_error(mpm_to_table(mat_u, matC = mat_c_na, xmax = xmax))
  expect_error(mpm_to_table(mat_u,
    matF = mat_f_named, start = "sm",
    xmax = xmax
  ))
})
