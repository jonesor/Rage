
test_that("vr functions work correctly", {
  
  library(popbio)
  
  ### test vr_U
  U1 <- rbind(c(  0,   0,   0,   0),
              c(0.5, 0.2, 0.1,   0),
              c(  0, 0.3, 0.3, 0.1),
              c(  0,   0, 0.5, 0.6))
  
  surv1 <- vr_survival(U1)
  expect_type(surv1, "double")
  expect_length(surv1, 1)
  
  surv2 <- vr_survival(U1, exclude_col = 4)
  expect_true(surv2 != surv1)
  
  grow1 <- vr_growth(U1)
  expect_type(grow1, "double")
  expect_length(grow1, 1)
  
  grow2 <- vr_growth(U1, posU = matrix(TRUE, nrow = 4, ncol = 4))
  expect_true(grow2 != grow1)
  
  shri1 <- vr_shrinkage(U1)
  expect_type(shri1, "double")
  expect_length(shri1, 1)
  
  shri2 <- vr_shrinkage(U1, exclude_row = 2:3)
  expect_true(is.na(shri2))
  
  stas1 <- vr_stasis(U1)
  expect_type(stas1, "double")
  expect_length(stas1, 1)
  
  stas2 <- vr_stasis(U1, posU = matrix(TRUE, nrow = 4, ncol = 4))
  expect_true(stas2 != stas1)
  
  dent1 <- vr_dorm_enter(U1, dorm_stages = 4)
  expect_type(dent1, "double")
  expect_length(dent1, 1)
  
  dexi1 <- vr_dorm_exit(U1, dorm_stages = 4)
  expect_type(dexi1, "double")
  expect_length(dexi1, 1)
  
  # reproduction
  R1 <- rbind(c(  0,   0, 1.1, 1.6),
              c(  0,   0, 0.8, 0.4),
              c(  0,   0,   0,   0),
              c(  0,   0,   0,   0))
  
  fecu1 <- vr_fecundity(U1, R1)
  expect_type(fecu1, "double")
  expect_length(fecu1, 1)
  
  # weighting
  v <- reproductive.value(U1 + R1)
  fecu2 <- vr_fecundity(U1, R1, weights_row = v)
  expect_true(fecu2 != fecu1)
  
  # no survival in a stage with reproduction
  U2 <- rbind(c(  0,   0,   0,   0),
              c(0.5, 0.2, 0.1,   0),
              c(  0, 0.3, 0.3,   0),
              c(  0,   0, 0.5,   0))
  
  fecu3 <- vr_fecundity(U2, R1)
  expect_true(fecu3 != fecu1)
})


test_that("vr functions warn and fail gracefully", {
  
  expect_error(vr_survival(mat_u_na))
  expect_error(vr_growth(mat_u_na))
  expect_error(vr_shrinkage(mat_u_na))
  expect_error(vr_stasis(mat_u_na))
  expect_error(vr_fecundity(mat_u_na))
})
