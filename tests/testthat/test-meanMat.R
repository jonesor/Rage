
test_that("meanMat works correctly", {
  s <- 3 # matrix dimension

  # generate lists of matrices
  x1 <- x2 <- replicate(5, matrix(runif(s^2), s, s), simplify = FALSE)
  x2[[1]][1, 1] <- NA

  y1 <- meanMat(x1)
  expect_true(inherits(y1, "matrix"))
  expect_true(nrow(y1) == s)

  # na handling
  y2 <- meanMat(x2)
  expect_true(inherits(y2, "matrix"))
  expect_true(is.na(y2[1, 1]))

  y3 <- meanMat(x2, na.rm = TRUE)
  expect_true(!is.na(y3[1, 1]))
})

test_that("meanMat warns and fails gracefully", {
  x3 <- list(matrix(0:3, nrow = 2), matrix(1:9, nrow = 3))
  expect_error(meanMat(x3))
})
