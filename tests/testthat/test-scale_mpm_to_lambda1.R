test_that("scale_mpm_to_lambda1 works correctly with matA only", {
  matA <- mat_u + mat_f

  x <- scale_mpm_to_lambda1(matA = matA)

  expect_type(x, "list")
  expect_length(x, 6L)
  expect_null(x$matU)
  expect_null(x$matF)
  expect_null(x$matC)
  expect_true(is.numeric(x$lambda_original))
  expect_equal(x$lambda_scaled, 1, tolerance = 1e-8)
  expect_equal(as.numeric(Re(popdemo::eigs(x$matA, what = "lambda"))),
    1,
    tolerance = 1e-8
  )
})

test_that("scale_mpm_to_lambda1 works correctly with component matrices", {
  matA <- mat_u + mat_f + mat_c

  x1 <- scale_mpm_to_lambda1(matU = mat_u, matF = mat_f)
  expect_type(x1, "list")
  expect_length(x1, 6L)
  expect_true(all(x1$matC == 0))
  expect_equal(x1$lambda_scaled, 1, tolerance = 1e-8)
  expect_equal(x1$matA, x1$matU + x1$matF + x1$matC, tolerance = 1e-8)

  x2 <- scale_mpm_to_lambda1(matU = mat_u, matC = mat_c)
  expect_true(all(x2$matF == 0))
  expect_equal(x2$lambda_scaled, 1, tolerance = 1e-8)
  expect_equal(x2$matA, x2$matU + x2$matF + x2$matC, tolerance = 1e-8)

  x3 <- scale_mpm_to_lambda1(matU = mat_u, matF = mat_f, matC = mat_c)
  expect_equal(x3$lambda_scaled, 1, tolerance = 1e-8)
  expect_equal(x3$matA, x3$matU + x3$matF + x3$matC, tolerance = 1e-8)

  x4 <- scale_mpm_to_lambda1(matU = mat_u, matF = mat_f, matC = mat_c, matA = matA)
  expect_equal(x4$lambda_scaled, 1, tolerance = 1e-8)
  expect_equal(x4$matA, x4$matU + x4$matF + x4$matC, tolerance = 1e-8)
})

test_that("scale_mpm_to_lambda1 preserves stage names", {
  matA_named <- mat_u_named + mat_f_named

  x <- scale_mpm_to_lambda1(matU = mat_u_named, matF = mat_f_named, matA = matA_named)

  expect_identical(dimnames(x$matA), dimnames(matA_named))
  expect_identical(dimnames(x$matU), dimnames(mat_u_named))
  expect_identical(dimnames(x$matF), dimnames(mat_f_named))
  expect_identical(dimnames(x$matC), dimnames(mat_u_named))
})

test_that("scale_mpm_to_lambda1 warns and fails gracefully", {
  expect_error(scale_mpm_to_lambda1())
  expect_error(scale_mpm_to_lambda1(matU = mat_u))
  expect_error(scale_mpm_to_lambda1(matF = mat_f))
  expect_error(scale_mpm_to_lambda1(matC = mat_c))
  expect_error(scale_mpm_to_lambda1(matU = mat_u, matF = mat_f_notsq))
  expect_error(scale_mpm_to_lambda1(matU = mat_u, matF = mat_f, matC = mat_c, matA = mat_u + mat_f))
  expect_error(scale_mpm_to_lambda1(matU = mat_u_na, matF = mat_f))
  expect_error(scale_mpm_to_lambda1(matA = matrix(0, nrow = 4, ncol = 4)))
  expect_error(scale_mpm_to_lambda1(matU = mat_u_named, matF = mat_f))
})
