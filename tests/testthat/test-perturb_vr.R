test_that("perturb_vr works correctly", {
  x <- perturb_vr(mat_u, mat_f)
  expect_type(x, "list")
  expect_equal(length(x), 5)
  expect_equal(x$clonality, 0)

  x_elast <- perturb_vr(mat_u, mat_f, type = "elasticity")
  expect_equal(x_elast$clonality, 0)

  x_f_zero <- suppressWarnings(
    perturb_vr(mat_u, matF = mat_f_zero, matC = mat_c)
  )
  expect_equal(x_f_zero$fecundity, 0)

  # check no growth/shrinkage in age-only model
  x_age <- perturb_vr(mat_u_age, mat_f_age)
  expect_equal(x_age$growth, 0)
  expect_equal(x_age$shrinkage, 0)

  # check works with custom demog_stat function
  fn_custom <- function(x) {
    return(1.0)
  }
  x_cust <- perturb_vr(mat_u, mat_f, demog_stat = "fn_custom")
  expect_true(all(x_cust == 0))
})

test_that("perturb_vr warns and fails gracefully", {
  notfn <- "rtpQswpAclFKhAmw"
  expect_error(perturb_vr(mat_u, mat_f, demog_stat = notfn))
})
