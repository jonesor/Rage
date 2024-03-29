test_that("entropy_k_stage works correctly with stage-based matrices", {
  data(leslie_mpm1)
  testMatrix <- mpm1$matU

  x <- entropy_k_stage(testMatrix)
  expect_true(inherits(x, "numeric"))
  expect_length(x, 1L)
  expect_gte(x, 0)

  x <- entropy_k_stage(testMatrix,
    init_distrib = runif(nrow(testMatrix))
  )
  expect_true(inherits(x, "numeric"))
  expect_length(x, 1L)
  expect_gte(x, 0)

  x <- entropy_k_stage(testMatrix, max_age = 5)
  expect_true(inherits(x, "numeric"))
  expect_length(x, 1L)
  expect_gte(x, 0)

  x <- entropy_k_stage(testMatrix, max_age = 500)
  expect_true(inherits(x, "numeric"))
  expect_length(x, 1L)
  expect_gte(x, 0)


  # Does function return NA when the matrix contains an NA?
  testMatrix_NA <- testMatrix
  testMatrix_NA[1, 1] <- NA
  x <- entropy_k_stage(testMatrix_NA)

  expect_true(is.na(x))
})

test_that("entropy_k_stage warns and fails gracefully with stage-based matrices", {
  data(leslie_mpm1)
  testMatrix <- mpm1$matU


  # Matrix not square
  testMatrix_ns <- testMatrix[-1]
  expect_error(entropy_k_stage(testMatrix_ns))

  # Matrix not numeric
  testMatrix_char <- apply(testMatrix, c(1, 2), as.character)
  expect_error(entropy_k_stage(testMatrix_char))
})
