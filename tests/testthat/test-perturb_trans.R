
test_that("perturb_trans works correctly", {
  
  x_sens <- perturb_trans(mat_u, mat_f)
  expect_type(x_sens, "list")
  expect_equal(length(x_sens), 5)
  expect_true(is.na(x_sens$clonality))
  
  x_elas <- perturb_trans(mat_u, mat_f, type = "elasticity")
  expect_type(x_elas, "list")
  expect_equal(length(x_elas), 5)
  expect_true(is.na(x_elas$clonality))
  
  x_f_zero <- suppressWarnings(
    perturb_trans(mat_u, matF = mat_f_zero, matC = mat_c)
  )
  expect_true(is.na(x_f_zero$fecundity))
  
  # check no stasis/retrogression in age-only model
  x_age <- perturb_trans(mat_u_age, mat_f_age)
  expect_true(is.na(x_age$stasis))
  expect_true(is.na(x_age$retro))
  
  # check works with custom demog_stat function
  fn_custom <- function(x) return(1)
  x_cust <- perturb_trans(mat_u, mat_f, mat_c, demog_stat = "fn_custom")
  expect_true(all(x_cust == 0))
})


test_that("perturb_trans warns and fails gracefully", {
  notfn <- "rtpsqwpaclfkhamw"
  expect_error(perturb_trans(mat_u, mat_f, demog_stat = notfn))
})
