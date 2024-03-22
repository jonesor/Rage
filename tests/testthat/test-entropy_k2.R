test_that("entropy_k2 works correctly", {
  data(leslie_mpm1)
  testMatrix <- leslie_mpm1$matU + leslie_mpm1$matF
  
  x <- entropy_k2(testMatrix)
  expect_true(inherits(x,"numeric"))
  expect_length(x, 1L)
  expect_gte(x, 0)
  
  #Does function return NA when the matrix contains an NA?
  testMatrix_NA <- testMatrix
  testMatrix_NA[1,1] <- NA
  x <- entropy_k2(testMatrix_NA)
  
  expect_true(is.na(x))
  
})

test_that("entropy_k2 warns and fails gracefully", {
  data(leslie_mpm1)
  testMatrix <- leslie_mpm1$matU + leslie_mpm1$matF
  
  
  # Matrix not square
  testMatrix_ns <- testMatrix[-1]
  expect_error(entropy_k2(testMatrix_ns))
  
  # Matrix not numeric
  testMatrix_char <- apply(testMatrix, c(1,2), as.character)
  expect_error(entropy_k2(testMatrix_char))
})
