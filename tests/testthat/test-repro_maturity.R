context("repro_maturity")

test_that("repro_maturity functions work correctly", {
  
  x_prob1 <- mature_prob(mat_u, mat_f)
  expect_is(x_prob1, "numeric")
  expect_length(x_prob1, 1L)
  
  x_prob2 <- mature_prob(mat_u, mat_f, startLife = 2)
  expect_gt(x_prob2, x_prob1)
  
  
  x_age1 <- mature_age(mat_u, mat_f)
  expect_is(x_age1, "numeric")
  expect_length(x_age1, 1L)
  
  x_age2 <- mature_age(mat_u, mat_f, startLife = 2)
  expect_lt(x_age2, x_age1)
  
  
  x_life1 <- mature_life_expect(mat_u, mat_f)
  expect_is(x_life1, "numeric")
  expect_length(x_life1, 1L)
})


test_that("repro_maturity functions warn and fail gracefully", {
  
  expect_warning(mature_prob(mat_u_survissue, mat_f))
  expect_error(mature_age(mat_u_na, mat_f))
  expect_error(mature_life_expect(mat_u, mat_f_na))
  expect_error(mature_life_expect(mat_u, mat_f, startLife = 10))
})

