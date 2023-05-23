test_that("mature_ functions work correctly", {
  x_prob1 <- mature_prob(mat_u, mat_f)
  expect_type(x_prob1, "double")
  expect_length(x_prob1, 1L)

  x_prob2 <- mature_prob(mat_u, mat_f, start = 2)
  expect_gt(x_prob2, x_prob1)

  x_age1 <- mature_age(mat_u, mat_f)
  expect_type(x_age1, "double")
  expect_length(x_age1, 1L)

  x_age2 <- mature_age(mat_u, mat_f, start = 2)
  expect_lt(x_age2, x_age1)

  # test life cycles with no connection from 'start' to repro stage(s)
  matU <- rbind(
    c(0.0, 0.0, 0.0),
    c(0.1, 0.2, 0.0),
    c(0.0, 0.0, 0.0)
  )

  matF <- rbind(
    c(0.0, 0.0, 5.3),
    c(0.0, 0.0, 0.0),
    c(0.0, 0.0, 0.0)
  )

  x_age3 <- mature_age(matU, matF, start = 1)
  expect_true(is.na(x_age3))
})


test_that("mature_ functions warn and fail gracefully", {
  expect_warning(mature_prob(mat_u_survissue, mat_f))
  expect_error(mature_age(mat_u_na, mat_f))
  expect_error(mature_prob(mat_u, mat_f, start = "stage name"))
  expect_error(mature_age(mat_u_named, mat_f_named, start = "invalid stage"))
})
