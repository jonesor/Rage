test_that("mature_distrib works correctly", {
  repstages <- repro_stages(mat_f)
  x1 <- mature_distrib(mat_u, repro_stages = repstages)
  expect_length(x1, ncol(mat_u))
  expect_identical(sum(x1), 1)
  expect_identical(x1[x1 > 0], 1) # repro maturity in only 1 stage class

  # test using named stages
  x2 <- mature_distrib(mat_u_named, repro_stages = c("lg", "xl"))
  expect_length(x2, ncol(mat_u_named))
  expect_identical(sum(x2), 1)
  expect_identical(x1, unname(x2))

  mat_u2 <- mat_u
  mat_u2[4, 2] <- 0.1
  x2 <- mature_distrib(mat_u2, repro_stages = repstages)
  expect_length(x2[x2 > 0], 2) # repro maturity in 2 stage classes

  # test life cycle with no connection from 'start' to repro stage(s)
  matU <- rbind(
    c(0.0, 0.0, 0.0),
    c(0.1, 0.2, 0.0),
    c(0.0, 0.0, 0.0)
  )

  x3 <- mature_distrib(matU, start = 1, repro_stages = c(FALSE, FALSE, TRUE))
  expect_true(all(x3 == 0))
})

test_that("mature_distrib warns and fails gracefully", {
  repstages <- repro_stages(mat_f)
  expect_error(mature_distrib(mat_u_na, repro_stages = repstages))
  expect_error(mature_distrib(mat_u, start = 10, repro_stages = repstages))
  expect_error(mature_distrib(mat_u_named, repro_stages = c(3, "lg")))
  expect_error(mature_distrib(mat_u_named, repro_stages = 3:6))
  expect_error(mature_distrib(mat_u_named, repro_stages = c("lg", "xl", "xxl")))
})
