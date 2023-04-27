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
  expect_identical(nrow(x_u), as.integer(xmax + 1))
  expect_identical(ncol(x_u), 7L)
  expect_identical(ncol(x_uf), 9L)
  expect_identical(ncol(x_uc), 9L)
  expect_identical(ncol(x_ufc), 13L)
  expect_identical(x_u, x_u_named)
  expect_identical(x_uf, x_uf_named)
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


mu_1 <- matrix(c(0.08,0.9,0.08,0.9),ncol = 2)
mu_1
mpm_to_table(mu_1)
expect_s3_class(mpm_to_table(mu_1), "data.frame")

expect_error(mpm_to_table(mu_1,remove_final = "Yes"))
expect_warning(mpm_to_table(mu_1,lx_crit = 0.1))

expect_s3_class(mpm_to_table(mu_1, radix = 1000), "data.frame")

expect_s3_class(mpm_to_table(mu_1, radix = 1000, 
                             remove_final = TRUE), "data.frame")

