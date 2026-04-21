test_that("mpm_collapse works correctly", {
  c1 <- list(1:2, 3:4)
  c2 <- list(1, 2, 3, 4) # should yield same as original
  c3 <- list(c("sm", "md"), c("lg", "xl")) # should yield same as c1
  c4 <- list(stage1 = c("sm", "md"), stage2 = c("lg", "xl"))

  x1 <- mpm_collapse(mat_u, mat_f, collapse = c1)
  x2 <- mpm_collapse(mat_u, mat_f, collapse = c2)
  x3 <- mpm_collapse(mat_u_named, mat_f_named, collapse = c3)
  x4 <- mpm_collapse(mat_u_named, mat_f_named, collapse = c4)

  expect_type(x1, "list")
  expect_length(x1, 4)
  expect_true(all(x1$matC == 0))

  expect_identical(ncol(x1$matA), length(c1))
  expect_identical(ncol(x2$matA), length(c2))
  expect_identical(ncol(x3$matA), length(c3))
  expect_identical(ncol(x4$matA), length(c4))

  expect_identical(x2$matU, mat_u)
  expect_identical(x2$matF, mat_f)

  expect_identical(x1, x3)

  expect_identical(colnames(x4$matA), names(c4))
})

test_that("mpm_collapse handles collapsed groups with zero SSD weights", {
  A <- matrix(c(
    0, 0.996209935189834, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0.991584639707654, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0.981367788179782, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0.959006021461268, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.395469318032286, 0, 0, 0, 0, 0.911050886198815, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0.427337118548192, 0, 0, 0, 0, 0,
    0.81275654994096, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.459469795638312, 0,
    0, 0, 0, 0, 0, 0.630395869496673, 0, 0, 0, 0, 0, 0, 0, 0,
    0.49155469074513, 0, 0, 0, 0, 0, 0, 0, 0.358122270025713, 0, 0, 0,
    0, 0, 0, 0, 0.523257244553273, 0, 0, 0, 0, 0, 0, 0, 0,
    0.101736892726593, 0, 0, 0, 0, 0, 0, 0.554226366194208, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0.00618162578413958, 0, 0, 0, 0, 0,
    0.584100587303554, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1.21341309051636e-05, 0, 0, 0, 0, 0.612514861948583, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1.1461992080595e-11, 0, 0, 0,
    0.639107841724659, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    4.47646892042388e-25, 0, 0, 0.663529428707613, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 6.45892956429893e-55, 0, 0.685448388953421, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.50204082860751e-121,
    0.704559797110107, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ), byrow = FALSE, ncol = 16)

  split_A <- mpm_split(A)
  stages <- list(1:2, 3:4, 5:6, 7:8, 9:10, 11:12, 13:14, 15:16)

  x <- mpm_collapse(
    matU = split_A$matU,
    matF = split_A$matF,
    collapse = stages
  )

  expect_false(any(is.nan(x$matU)))
  expect_false(any(is.nan(x$matA)))
  expect_equal(x$matU[7:8, 8], c(2.23823446021194e-25, 1.38003227351482e-79))
})

test_that("mpm_collapse warns and fails gracefully", {
  expect_error(mpm_collapse(mat_u_na, mat_f, collapse = list(1:2, 3:4)))
  expect_error(mpm_collapse(mat_u_named, mat_f_named, collapse = list(1:4, 5)))
  expect_error(mpm_collapse(mat_u_named, mat_f_named,
    collapse = list(c("sm", "md", "lg", "xl"), "xxl")
  ))
  expect_error(mpm_collapse(mat_u_named, mat_f_named,
    collapse = list(1:3, "xl")
  ))
})
