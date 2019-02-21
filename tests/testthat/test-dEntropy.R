context("dEntropy")

test_that("dEntropy works correctly", {
  
  x <- dEntropy(mat_u, mat_f)

  expect_length(x, 1L)
  expect_true(x > 0)
})

test_that("dEntropy warns and fails gracefully", {
  
  expect_warning(dEntropy(mat_u_survissue, mat_f))
  expect_error(dEntropy(mat_u_na, mat_f))
  expect_error(dEntropy(mat_u, mat_f_na))
  expect_error(dEntropy(mat_u, mat_f, startLife = 10))
})
