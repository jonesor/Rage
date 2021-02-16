
test_that("mpm_split works correctly", {
  
  matA <- mat_u + mat_f
  
  x <- mpm_split(matA)
  
  expect_type(x, "list")
  expect_length(x, 3)
  expect_true(all(x$matC == 0))
  
  expect_equal(ncol(mat_u), ncol(x$matU))
})
