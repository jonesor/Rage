context("checkValidMat")

test_that("checkValidMat works correctly", {
  
  matU <- matrix(1:4, nrow = 2)
  expect_silent(checkValidMat(matU))
  
  # not matrix
  matU <- list(1:4)
  expect_error(checkValidMat(M = matU))
  
  # not numeric
  matU <- matrix(letters[1:4], nrow = 2)
  expect_error(checkValidMat(M = matU))
  
  # not square
  matU <- matrix(rep(0, 6), nrow = 2)
  expect_error(checkValidMat(M = matU))
  
  # contains all NA
  matU <- matrix(NA_real_, nrow = 2, ncol = 2)
  expect_error(checkValidMat(matU))
  
  # contains NA
  matU <- matrix(c(1, NA, 3, 4), nrow = 2)
  expect_error(checkValidMat(matU))
  
  # survival issue
  matU <- matrix(c(1, 1, 1, 1), nrow = 2)
  expect_warning(checkValidMat(matU, warn_surv_issue = TRUE))
  
  # all zeros
  matU <- matrix(rep(0, 4), nrow = 2)
  expect_warning(checkValidMat(M = matU))
})

