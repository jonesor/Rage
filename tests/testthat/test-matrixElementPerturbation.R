context("matrixElementPerturbation")

test_that("matrixElementPerturbation works correctly", {
  
  x <- matrixElementPerturbation(mat_u, mat_f)
  expect_is(x, "data.frame")
  expect_equal(ncol(x), 10)
  expect_equal(x$SClonality, 0)
  expect_equal(x$EClonality, 0)
  
  x_f_zero <- suppressWarnings(
    matrixElementPerturbation(mat_u, matF = mat_f_zero, matC = mat_c)
  )
  expect_equal(x_f_zero$SFecundity, 0)
  expect_equal(x_f_zero$EFecundity, 0)
  
  # check no stasis/retrogression in age-only model
  x_age <- matrixElementPerturbation(mat_u_age, mat_f_age)
  expect_equal(x_age$SStasis, 0)
  expect_equal(x_age$SRetrogression, 0)
  expect_equal(x_age$EStasis, 0)
  expect_equal(x_age$ERetrogression, 0)
  
  # check works with custom demogstat function
  fn_custom <- function(x) return(1)
  x_cust <- matrixElementPerturbation(mat_u, mat_f, demogstat = "fn_custom")
  expect_true(all(x_cust == 0))
})

test_that("matrixElementPerturbation warns and fails gracefully", {
  notfn <- "rtpsqwpaclfkhamw"
  expect_error(matrixElementPerturbation(mat_u, mat_f, demogstat = notfn))
})
