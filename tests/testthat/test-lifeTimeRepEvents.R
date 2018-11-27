context("lifeTimeRepEvents")

test_that("lifeTimeRepEvents works correctly", {
  
  x <- lifeTimeRepEvents(mat_u, mat_f)
  expect_is(x, "list")
  expect_length(x, 4L)
  
  y <- lifeTimeRepEvents(mat_u, mat_f, startLife = 2)
  expect_is(y, "list")
  expect_length(y, 4L)
  expect_gt(y$remainingMatureLifeExpectancy, x$remainingMatureLifeExpectancy)
})

test_that("lifeTimeRepEvents warns and fails gracefully", {
  
  expect_warning(lifeTimeRepEvents(mat_u_survissue, mat_f))
  expect_error(lifeTimeRepEvents(mat_u_na, mat_f))
  expect_error(lifeTimeRepEvents(mat_u, mat_f_na))
  expect_error(lifeTimeRepEvents(mat_u, mat_f, startLife = 10))
})
