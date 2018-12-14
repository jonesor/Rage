context("splitMatrix")

test_that("splitMatrix works correctly", {
  
  matA <- mat_u + mat_f
  
  x <- splitMatrix(matA)
  
  expect_is(x, "list")
  expect_length(x, 3)
  expect_true(all(x$matC == 0))
  
  expect_equal(ncol(mat_u), ncol(x$matU))
})
