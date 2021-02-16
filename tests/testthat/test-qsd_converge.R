
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

test_that("qsd_converge works w/ non-ergodic matrix", {
  
  mat_no_ergo <- matrix(
    c(
      0,   2.2,  5.8,  0,    0,
      0.8, 0,    0,    0,    0,
      0,   0.89, 0,    0,    0,
      0,   0,    0.93, 0,    0,
      0,   0,    0,    0.75, 0.5
    ),
    nrow = 5, ncol = 5, byrow = TRUE
  )
  
  # fast convergence
  t_qsd <- qsd_converge(mat_no_ergo)
  
  expect_length(t_qsd, 1L)
  expect_equal(t_qsd, 10L)
  
  # multi-state w/ non-ergodic also works.
  f_qsd <- qsd_converge(mat_no_ergo,
                        start = c(1, 2, 3, 2, 0))
  
  expect_length(f_qsd, 1L)
  expect_equal(f_qsd, 8L)
  
})
