test_that("shape_rep works correctly", {
  # constant mx
  mx1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1)
  s1a <- shape_rep(mx1)
  expect_equal(s1a, 0)

  # constant hazard, custom xlim
  s1b <- shape_rep(mx1, xmin = 2, xmax = 5)
  expect_equal(s1b, 0)

  # decreasing mx
  mx2 <- seq(1, 0, -0.1)
  s2 <- shape_rep(mx2)
  expect_true(s2 > 0 & s2 < 0.5)

  # increasing mx
  mx3 <- seq(0, 1, 0.1)
  s3 <- shape_rep(mx3)
  expect_true(s3 < 0 & s3 > -0.5)

  # check works with data frame
  lt <- data.frame(x = seq_along(mx3) - 1, mx = mx3)
  s4 <- shape_rep(lt)
  expect_equal(s4, s3)
})


test_that("shape_rep warns and fails gracefully", {
  # negative reproduction
  # You appear to have minus-babies (check mx)
  expect_error(shape_rep(c(0, 1, 1, 1, -0.1, 0)))

  # < 3 nozero values of mx
  # must have > 2 nonzero values of mx to calculate shape
  expect_error(shape_rep(c(0, 0.5, 0.5)))

  # "'rep' doesn't contain both x and mx"
  rep1 <- list(mx = c(0, 0, 0.3, 0.4, 0.5, 0.6))
  expect_error(shape_rep(rep1))

  rep1 <- list(years = 0:5, mx = c(0, 0, 0.3, 0.4, 0.5, 0.6))
  expect_error(shape_rep(rep1))

  # x and mx must be the same length
  rep1 <- list(x = 0:4, mx = c(0, 0, 0.3, 0.4, 0.5, 0.6))
  expect_error(shape_rep(rep1))

  # much as we'd like to reverse ageing, x must all be ascending
  rep1 <- list(x = c(0, 1, 2, 3, 4, 3), mx = c(0, 0, 0.3, 0.4, 0.5, 0.6))
  expect_error(shape_rep(rep1))
})
