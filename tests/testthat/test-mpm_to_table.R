context("mpm_to_table")

test_that("mpm_to_table works correctly", {
  
  xmax <- 20
  
  x_u <- mpm_to_table(mat_u, xmax = xmax, lxCrit = 0)
  x_uf <- mpm_to_table(mat_u, mat_f, xmax = xmax, lxCrit = 0)
  x_uc <- mpm_to_table(mat_u, matC = mat_c, xmax = xmax, lxCrit = 0)
  x_ufc <- mpm_to_table(mat_u, mat_f, mat_c, xmax = xmax, lxCrit = 0)
  
  expect_is(x_u, "data.frame")
  expect_equal(nrow(x_u), xmax + 1)
  expect_equal(ncol(x_u), 7)
  expect_equal(ncol(x_uf), 9)
  expect_equal(ncol(x_uc), 9)
  expect_equal(ncol(x_ufc), 13)
})

test_that("mpm_to_table warns and fails gracefully", {
  
  xmax <- 20
  
  expect_error(mpm_to_table(mat_u_na, mat_f, xmax = xmax))
  expect_error(mpm_to_table(mat_u, mat_f_na, xmax = xmax))
  expect_error(mpm_to_table(mat_u, matC = mat_c_na, xmax = xmax))
  expect_warning(mpm_to_table(mat_u_survissue, xmax = xmax))
  expect_warning(mpm_to_table(mat_u, mat_f_zero, xmax = xmax))
  expect_warning(mpm_to_table(mat_u, matC = mat_c_zero, xmax = xmax))
})
