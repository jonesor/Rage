context("shape_surv")

test_that("shape_surv works correctly", {
  
  # constant hazard
  lx1 <- 0.5^(0:20)
  s1a <- shape_surv(lx1)
  expect_equal(s1a, 0)
  
  # constant hazard, custom xlim
  s1b <- shape_surv(lx1, xmin = 2, xmax = 10)
  expect_equal(s1b, 0)
  
  # increasing hazard
  x2 <- seq(0, 1, 0.1)
  hx2 <- 0.5 * x2
  lx2 <- hx_to_lx(hx2)
  s2 <- shape_surv(lx2)
  expect_true(s2 > 0 & s2 < 0.5)
  
  # declining hazard
  x3 <- seq(0, 1, 0.1)
  hx3 <- 0.2 - 0.2 * x3
  lx3 <- hx_to_lx(hx3)
  s3 <- shape_surv(lx3)
  expect_true(s3 < 0 & s3 > -0.5)
  
  # check works with data frame
  lt <- data.frame(x = x3, lx = lx3)
  s4 <- shape_surv(lt)
  expect_equal(s4, s3)
})


test_that("shape_surv warns and fails gracefully", {
  
  # zombies
  expect_error(shape_surv(c(1, 0.7, 0.8, 0.3)))
  
  # < 3 nozero values of lx
  expect_error(shape_surv(c(1, 0.5, 0)))
  
  # lx[1] != 0
  expect_error(shape_surv(c(0.8, 0.7, 0.6)))
})