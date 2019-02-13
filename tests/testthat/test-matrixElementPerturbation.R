context("perturb_trans")

test_that("perturb_trans works correctly", {
  
  x <- perturb_trans(mat_u, mat_f)
  expect_is(x, "list")
  expect_equal(length(x), 5)
  expect_true(is.na(x$clonality))
  
  x_f_zero <- suppressWarnings(
    perturb_trans(mat_u, matF = mat_f_zero, matC = mat_c)
  )
  expect_true(is.na(x_f_zero$fecundity))
  
  # check no stasis/retrogression in age-only model
  x_age <- perturb_trans(mat_u_age, mat_f_age)
  expect_true(is.na(x_age$stasis))
  expect_true(is.na(x_age$retro))
  
  # check works with custom demogstat function
  fn_custom <- function(x) return(1)
  x_cust <- perturb_trans(mat_u, mat_f, mat_c, demogstat = "fn_custom")
  expect_true(all(x_cust == 0))
})


test_that("perturb_trans warns and fails gracefully", {
  notfn <- "rtpsqwpaclfkhamw"
  expect_error(perturb_trans(mat_u, mat_f, demogstat = notfn))
})
