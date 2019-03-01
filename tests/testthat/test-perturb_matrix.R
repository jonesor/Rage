context("perturb_matrix")

test_that("perturb_matrix works correctly", {
  
  mat_a <- mat_u + mat_f
  
  x_sens1 <- perturb_matrix(mat_a)
  expect_is(x_sens1, "matrix")
  expect_true(nrow(x_sens1) == nrow(mat_a) & ncol(x_sens1) == ncol(mat_a))
  
  x_sens2 <- perturb_matrix(mat_a, pert = 0.2)
  expect_is(x_sens1, "matrix")
  expect_true(sum(x_sens1) != sum(x_sens2))
  
  x_elas1 <- perturb_matrix(mat_a, type = "elasticity")
  expect_is(x_elas1, "matrix")
  expect_true(abs(sum(x_elas1) - 1) < 1e-6)
  
  # check works with custom demogstat function
  fn_custom <- function(x) return(1)
  x_cust <- perturb_matrix(mat_a, demogstat = "fn_custom")
  expect_true(all(x_cust == 0))
})


test_that("perturb_matrix warns and fails gracefully", {
  notfn <- "rtpsqwpaclfkhamw"
  expect_error(perturb_matrix(mat_a, demogstat = notfn))
})

