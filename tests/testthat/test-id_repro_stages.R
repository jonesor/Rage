context("id_repro_stages")

test_that("id_repro_stages works correctly", {

  x <- id_repro_stages(mat_f)
  
  expect_is(x, "logical")
  expect_length(x, ncol(mat_f))
  
  x_na_1 <- id_repro_stages(mat_f_na, na.handling = "return.true")
  x_na_2 <- id_repro_stages(mat_f_na, na.handling = "return.na")
  x_na_3 <- id_repro_stages(mat_f_na, na.handling = "return.false")
  
  expect_false(any(is.na(x_na_1)))
  expect_true(any(is.na(x_na_2)))
  expect_false(any(is.na(x_na_3)))
})

test_that("id_repro_stages warns and fails gracefully", {
  
  m_na <- matrix(NA_real_, nrow = 4, ncol = 4) 
  
  expect_error(id_repro_stages(m_na))
  expect_error(id_repro_stages(mat_f, na.handling = "eigen"))
})
