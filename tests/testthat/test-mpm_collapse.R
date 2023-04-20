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

  expect_equal(ncol(x1$matA), length(c1))
  expect_equal(ncol(x2$matA), length(c2))
  expect_equal(ncol(x3$matA), length(c3))
  expect_equal(ncol(x4$matA), length(c4))

  expect_equal(x2$matU, mat_u)
  expect_equal(x2$matF, mat_f)

  expect_identical(x1, x3)

  expect_equal(colnames(x4$matA), names(c4))
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
