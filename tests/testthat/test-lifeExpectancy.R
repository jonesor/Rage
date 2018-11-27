context("lifeExpectancy")

test_that("lifeExpectancy works correctly", {
  
  l0 <- lifeExpectancy(mat_u)
  l0_singular <- lifeExpectancy(mat_u_singular)

  expect_length(l0, 1L)
  expect_true(l0 > 0)
  expect_equal(l0_singular, NA_real_)
})

test_that("lifeExpectancy warns and fails gracefully", {
  
  expect_warning(lifeExpectancy(mat_u_zero))
  expect_warning(lifeExpectancy(mat_u_survissue))
  expect_error(lifeExpectancy(mat_u, startLife = 10))
  expect_error(lifeExpectancy(mat_u_na))
})
