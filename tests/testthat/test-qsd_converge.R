context("qsd_converge")

test_that("qsd_converge works correctly", {
  
  x1 <- qsd_converge(mat_u)
  expect_length(x1, 1L)
  expect_true(x1 > 0)
  
  # matU all zero
  x2 <- suppressWarnings(qsd_converge(mat_u_zero))
  expect_equal(x2, 1L)
  
  # multiple starting stages
  x3 <- qsd_converge(mat_u, start = c(0.8, 0.2, 0, 0))
  expect_length(x3, 1L)
  expect_true(x3 > 0)
})

test_that("qsd_converge warns and fails gracefully", {
  expect_error(qsd_converge(mat_u_na))
  expect_error(qsd_converge(mat_u, start = 10))
})
