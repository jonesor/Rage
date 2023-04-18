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

test_that("checkValidMat produces warnings in mpm_to_table", {
  xmax <- 20
  surv_issue_regex <- "^Argument mat_u_survissue has at least one stage-specific survival*"
  repro_regex <- "All elements of mat_f_zero are zero"
  clone_regex <- "All elements of mat_c_zero are zero"

  expect_warning(checkValidMat(mat_u_survissue, warn_surv_issue = TRUE),
    regexp = surv_issue_regex
  )
  expect_warning(checkValidMat(mat_f_zero),
    regexp = repro_regex
  )
  expect_warning(checkValidMat(mat_c_zero),
    regexp = clone_regex
  )
})
