context("vr_vec")

test_that("vr_vec functions work correctly", {
  
  library(popbio)
  
  ### test vr_vec_U
  U1 <- rbind(c(  0,   0,   0,   0),
              c(0.5, 0.2, 0.1,   0),
              c(  0, 0.3, 0.3, 0.1),
              c(  0,   0, 0.5, 0.6))
  
  surv1 <- vr_vec_survival(U1)
  expect_is(surv1, "numeric")
  expect_length(surv1, 4)
  
  surv2 <- vr_vec_survival(U1, exclude_col = 4)
  expect_true(is.na(surv2[4]))
  
  surv3 <- vr_vec_survival(U1, exclude_col = c(FALSE, FALSE, FALSE, TRUE))
  expect_true(is.na(surv3[4]))
  
  grow1 <- vr_vec_growth(U1)
  expect_is(grow1, "numeric")
  expect_length(grow1, 4)
  
  grow2 <- vr_vec_growth(U1, posU = matrix(TRUE, nrow = 4, ncol = 4))
  expect_true(!is.na(grow2[1]))
  
  grow3 <- vr_vec_growth(U1, surv_only_na = FALSE)
  expect_true(!is.na(grow3[1]))
  
  shri1 <- vr_vec_shrinkage(U1)
  expect_is(shri1, "numeric")
  expect_length(shri1, 4)
  
  shri2 <- vr_vec_shrinkage(U1, exclude_row = 2)
  expect_true(is.na(shri2[3]))
  
  stas1 <- vr_vec_stasis(U1)
  expect_is(stas1, "numeric")
  expect_length(stas1, 4)
  
  stas2 <- vr_vec_stasis(U1, posU = matrix(TRUE, nrow = 4, ncol = 4))
  expect_equal(stas2[1], 0)
  
  dent1 <- vr_vec_dorm_enter(U1, dorm_stages = 4)
  expect_is(dent1, "numeric")
  expect_length(dent1, 4)
  expect_true(!is.na(dent1[3]))
  
  dexi1 <- vr_vec_dorm_exit(U1, dorm_stages = 4)
  expect_is(dexi1, "numeric")
  expect_length(dexi1, 4)
  expect_true(!is.na(dexi1[4]))
  
  
  R1 <- rbind(c(  0,   0, 1.1, 1.6),
              c(  0,   0, 0.8, 0.4),
              c(  0,   0,   0,   0),
              c(  0,   0,   0,   0))
  
  fecu1 <- vr_vec_reproduction(U1, R1)
  expect_is(fecu1, "numeric")
  expect_length(fecu1, 4)
  expect_true(all(is.na(fecu1[1:2])))
  
  # weighting
  v <- reproductive.value(U1 + R1)
  fecu2 <- vr_vec_reproduction(U1, R1, weights_row = v)
  expect_true(fecu2[3] != fecu1[3])
  
  # no survival in a stage with reproduction
  U2 <- rbind(c(  0,   0,   0,   0),
              c(0.5, 0.2, 0.1,   0),
              c(  0, 0.3, 0.3,   0),
              c(  0,   0, 0.5,   0))
  
  fecu3 <- vr_vec_reproduction(U2, R1)
  expect_equal(fecu3[4], sum(R1[,4]))
})


test_that("vr_vec functions warn and fail gracefully", {
  
  expect_error(vr_vec_survival(mat_u_na))
  expect_error(vr_vec_growth(mat_u_na))
  expect_error(vr_vec_shrinkage(mat_u_na))
  expect_error(vr_vec_stasis(mat_u_na))
  expect_error(vr_vec_reproduction(mat_u_na))
})
