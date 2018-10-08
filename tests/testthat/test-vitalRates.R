context("vitalRates")

test_that("vitalRates works correctly", {
  
  x <- vitalRates(mat_u, mat_f, mat_c, splitStages = "all")
  expect_is(x, "list")
  expect_length(x, 5)
  
  x <- vitalRates(mat_u, mat_f, mat_c, splitStages = "all", weights = "SSD")
  expect_length(x, 5)
  
  x <- vitalRates(mat_u, mat_f, splitStages = "ontogeny")
  expect_length(x, 9)
  
  ms <- c("prop", "active", "active", "active")
  x <- vitalRates(mat_u, mat_f, splitStages = "matrixStages", matrixStages = ms)
  expect_length(x, 10)
})

test_that("vitalRates warns and fails gracefully", {
  
  # weights has incorrect dimension
  expect_error(vitalRates(mat_u, mat_f, weights = c(0.5, 0.5)))
  
  # invalid splitStages
  expect_error(vitalRates(mat_u, mat_f, splitStages = "prop"))
  
  # matrixStages not provided when splitStages = "matrixStages"
  expect_error(vitalRates(mat_u, mat_f, splitStages = "matrixStages"))
  
  # matrixStages has incorrect dimension
  expect_error(vitalRates(mat_u, mat_f, splitStages = "matrixStages",
                          matrixStages = c("prop", "active")))
})
