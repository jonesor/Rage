
test_that("life_expect works correctly", {
  
  l0 <- life_expect(mat_u)
  l0_singular <- life_expect(mat_u_singular)

  expect_length(l0, 1L)
  expect_true(l0 > 0)
  expect_equal(l0_singular, NA_real_)
})

test_that("life_expect warns and fails gracefully", {
  
  expect_warning(life_expect(mat_u_zero))
  expect_warning(life_expect(mat_u_survissue))
  expect_error(life_expect(mat_u, start = 10))
  expect_error(life_expect(mat_u_na))
})
