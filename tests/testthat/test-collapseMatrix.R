context("mpm_collapse")

test_that("mpm_collapse works correctly", {
  
  c1 <- list(1:2, 3:4)
  c2 <- list(1, 2, 3, 4) # should yield same as original
  
  x1 <- mpm_collapse(mat_u, mat_f, collapse = c1)
  x2 <- mpm_collapse(mat_u, mat_f, collapse = c2)
  
  expect_is(x1, "list")
  expect_length(x1, 4)
  expect_true(all(x1$matC == 0))
  
  expect_equal(ncol(x1$matA), length(c1))
  expect_equal(ncol(x2$matA), length(c2))
  
  expect_equal(x2$matU, mat_u)
  expect_equal(x2$matF, mat_f)
})

test_that("mpm_collapse warns and fails gracefully", {
  
  expect_error(mpm_collapse(mat_u_na, mat_f, collapse = list(1:2, 3:4)))
})
