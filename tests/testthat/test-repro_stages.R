test_that("repro_stages works correctly", {
  x <- repro_stages(mat_f)

  expect_type(x, "logical")
  expect_length(x, ncol(mat_f))

  x_na_1 <- repro_stages(mat_f_na, na_handling = "return.true")
  x_na_2 <- repro_stages(mat_f_na, na_handling = "return.na")
  x_na_3 <- repro_stages(mat_f_na, na_handling = "return.false")

  expect_false(anyNA(x_na_1))
  expect_true(anyNA(x_na_2))
  expect_false(anyNA(x_na_3))
})

test_that("repro_stages warns and fails gracefully", {
  m_na <- matrix(NA_real_, nrow = 4, ncol = 4)

  expect_error(repro_stages(m_na))
  expect_error(repro_stages(mat_f, na_handling = "eigen"))
})
