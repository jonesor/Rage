context("makeLifeTable")

test_that("makeLifeTable works correctly", {
  
  nSteps <- 20
  
  x_u <- makeLifeTable(mat_u, nSteps = nSteps)
  x_uf <- makeLifeTable(mat_u, mat_f, nSteps = nSteps)
  x_uc <- makeLifeTable(mat_u, matC = mat_c, nSteps = nSteps)
  x_ufc <- makeLifeTable(mat_u, mat_f, mat_c, nSteps = nSteps)
  
  expect_is(x_u, "data.frame")
  expect_equal(nrow(x_u), nSteps)
  expect_equal(ncol(x_u), 7)
  expect_equal(ncol(x_uf), 9)
  expect_equal(ncol(x_uc), 9)
  expect_equal(ncol(x_ufc), 13)
})

test_that("makeLifeTable warns and fails gracefully", {
  
  nSteps <- 20
  
  expect_error(makeLifeTable(mat_u_na, mat_f, nSteps = nSteps))
  expect_error(makeLifeTable(mat_u, mat_f_na, nSteps = nSteps))
  expect_error(makeLifeTable(mat_u, matC = mat_c_na, nSteps = nSteps))
  expect_warning(makeLifeTable(mat_u_survissue, nSteps = nSteps))
  expect_warning(makeLifeTable(mat_u, mat_f_zero, nSteps = nSteps))
  expect_warning(makeLifeTable(mat_u, matC = mat_c_zero, nSteps = nSteps))
})
