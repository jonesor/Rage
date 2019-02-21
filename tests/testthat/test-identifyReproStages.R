context("identifyReproStages")

test_that("identifyReproStages works correctly", {

  x <- identifyReproStages(mat_f)
  
  expect_is(x, "logical")
  expect_length(x, ncol(mat_f))
  
  x_na_1 <- identifyReproStages(mat_f_na, na.handling = "return.true")
  x_na_2 <- identifyReproStages(mat_f_na, na.handling = "return.na")
  x_na_3 <- identifyReproStages(mat_f_na, na.handling = "return.false")
  
  expect_false(any(is.na(x_na_1)))
  expect_true(any(is.na(x_na_2)))
  expect_false(any(is.na(x_na_3)))
})

test_that("identifyReproStages warns and fails gracefully", {
  
  m_na <- matrix(NA_real_, nrow = 4, ncol = 4) 
  
  expect_error(identifyReproStages(m_na))
  expect_error(identifyReproStages(mat_f, na.handling = "eigen"))
})
