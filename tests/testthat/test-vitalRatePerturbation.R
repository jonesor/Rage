context("vitalRatePerturbation")


test_that("vitalRatePerturbation works correctly", {
  
  x <- vitalRatePerturbation(mat_u, mat_f)
  expect_is(x, "data.frame")
  expect_equal(ncol(x), 10)
  expect_equal(x$SClonality, 0)
  expect_equal(x$EClonality, 0)
  
  x_f_zero <- suppressWarnings(
    vitalRatePerturbation(mat_u, matF = mat_f_zero, matC = mat_c)
  )
  expect_equal(x_f_zero$SReproduction, 0)
  expect_equal(x_f_zero$EReproduction, 0)
  
  # check no growth/shrinkage in age-only model
  x_age <- vitalRatePerturbation(mat_u_age, mat_f_age)
  expect_equal(x_age$SGrowth, 0)
  expect_equal(x_age$SShrinkage, 0)
  expect_equal(x_age$EGrowth, 0)
  expect_equal(x_age$EShrinkage, 0)
  
  # check works with custom demogstat function
  fn_custom <- function(x) return(1.0)
  x_cust <- vitalRatePerturbation(mat_u, mat_f, demogstat = "fn_custom")
  expect_true(all(x_cust == 0))
})

test_that("vitalRatePerturbation warns and fails gracefully", {
  notfn <- "rtpQswpAclFKhAmw"
  expect_error(vitalRatePerturbation(mat_u, mat_f, demogstat = notfn))
})
