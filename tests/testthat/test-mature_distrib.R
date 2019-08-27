context("mature_distrib")

test_that("mature_distrib works correctly", {
  
  repstages <- repro_stages(mat_f)
  x1 <- mature_distrib(mat_u, repro_stages = repstages)
  expect_length(x1, ncol(mat_u))
  expect_equal(sum(x1), 1)
  expect_length(x1[x1 > 0], 1) # repro maturity in only 1 stage class
  
  mat_u2 <- mat_u
  mat_u2[4,2] <- 0.1
  x2 <- mature_distrib(mat_u2, repro_stages = repstages)
  expect_length(x2[x2 > 0], 2) # repro maturity in 2 stage classes
})

test_that("mature_distrib warns and fails gracefully", {
  
  repstages <- repro_stages(mat_f)
  expect_error(mature_distrib(mat_u_na, repro_stages = repstages))
  expect_error(mature_distrib(mat_u, start = 10, repro_stages = repstages))
})
