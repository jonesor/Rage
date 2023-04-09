test_that("vr_mat functions work correctly", {
  ### test vr_mat_U
  U1 <- rbind(
    c(0, 0, 0, 0),
    c(0.5, 0.2, 0.1, 0),
    c(0, 0.3, 0.3, 0.1),
    c(0, 0, 0.5, 0.6)
  )

  x1 <- vr_mat_U(U1)
  expect_true(inherits(U1, "matrix"))
  expect_true(all(is.na(x1[, 1])))

  x2 <- vr_mat_U(U1, surv_only_na = FALSE)
  expect_true(!all(is.na(x2[, 1])))

  # all transitions possible
  x3 <- vr_mat_U(U1, posU = matrix(TRUE, nrow = 4, ncol = 4))
  expect_true(!any(is.na(x3)))


  ### test vr_mat_R
  R1 <- rbind(
    c(0, 0, 1.1, 1.6),
    c(0, 0, 0.8, 0.4),
    c(0, 0, 0, 0),
    c(0, 0, 0, 0)
  )

  x4 <- vr_mat_R(U1, R1)
  expect_true(inherits(x4, "matrix"))
  expect_true(all(is.na(x4[, 1:2])))

  # all transitions possible
  x5 <- vr_mat_R(U1, R1, posR = matrix(TRUE, nrow = 4, ncol = 4))
  expect_true(!any(is.na(x5)))
})


test_that("vr_mat functions warn and fail gracefully", {
  expect_error(vr_mat_U(mat_u_na))
  expect_error(vr_mat_R(mat_u, mat_f_na))
})
