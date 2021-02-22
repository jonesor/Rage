
test_that("vitalRates functions work correctly", {
  
  ### test vitalRates
vrOut <- vitalRates(matU = mat_u,matF = mat_f)

expect_type(vrOut, "list")
expect_length(vrOut, 5)
expect_equal(vrOut$clo, 0)

vrOut_c <- vitalRates(matU = mat_u,matF = mat_f,matC = mat_c)

expect_type(vrOut_c, "list")
expect_length(vrOut_c, 5)


})


test_that("vitalRates functions warn and fail gracefully", {
  
  expect_error(vitalRates(mat_u_na, mat_f_na))
  expect_error(vitalRates(mat_u, mat_f_na))
  expect_error(vitalRates(mat_u_na, mat_f))
})
