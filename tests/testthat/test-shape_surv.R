test_that("shape_surv works correctly", {
  # constant hazard
  lx1 <- 0.5^(0:20)
  s1a <- shape_surv(lx1, trunc = TRUE)
  expect_identical(s1a, 0)

  # constant hazard, custom xlim
  s1b <- shape_surv(lx1, xmin = 2, xmax = 10, trunc = TRUE)
  expect_identical(s1b, 0)

  # increasing hazard
  x2 <- seq(0, 1, 0.1)
  hx2 <- 0.5 * x2
  lx2 <- hx_to_lx(hx2)
  s2 <- shape_surv(lx2, trunc = TRUE)
  expect_true(s2 > 0 & s2 < 0.5)

  # declining hazard
  x3 <- seq(0, 1, 0.1)
  hx3 <- 0.2 - 0.2 * x3
  lx3 <- hx_to_lx(hx3)
  s3 <- shape_surv(lx3, trunc = TRUE)
  expect_true(s3 < 0 & s3 > -0.5)

  # check works with data frame
  lt <- data.frame(x = x3, lx = lx3)
  s4 <- shape_surv(lt, trunc = TRUE)
  expect_identical(s4, s3)
})


test_that("shape_surv warns and fails gracefully", {
  # first element of lx is not 1
  expect_error(shape_surv(c(0.8, 0.7, 0.6), trunc = TRUE))

  # Zero not dealt with
  expect_error(shape_surv(c(1, 0.5, 0.25, 0)))

  # zombies
  expect_error(shape_surv(c(1, 0.7, 0.8, 0.3), trunc = TRUE))

  # < 3 nozero values of lx
  expect_error(shape_surv(c(1, 0.5, 0), trunc = TRUE))

  #' surv' doesn't contain both x and lx
  surv1 <- list(years = 0:3, lx = c(1, 0.8, 0.7, 0.6))
  expect_error(shape_surv(surv1))

  # x and lx must be the same length
  surv2 <- list(x = 0:2, lx = c(1, 0.8, 0.7, 0.6))
  expect_error(shape_surv(surv2))

  # lx must start with 1 where x[1] is 0
  surv3 <- list(x = 0:3, lx = c(0.9, 0.8, 0.7, 0.6))
  expect_error(shape_surv(surv3))

  # much as we'd like to reverse ageing, x must all be ascending
  surv3 <- list(x = c(0, 1, 2, 1), lx = c(1, 0.8, 0.7, 0.6))
  expect_error(shape_surv(surv3))
})
