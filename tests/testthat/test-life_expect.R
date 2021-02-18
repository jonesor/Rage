
test_that("life_expect works correctly", {
  
  l0 <- life_expect(mat_u)
  l0_singular <- life_expect(mat_u_singular)
  l1 <- life_expect(mat_u_named, start = "sm")

  expect_length(l0, 1L)
  expect_true(l0 > 0)
  expect_equal(l0_singular, NA_real_)
  expect_equal(l0, l1)
})

test_that("life_expect warns and fails gracefully", {
  
  expect_warning(life_expect(mat_u_zero))
  expect_warning(life_expect(mat_u_survissue))
  expect_error(life_expect(mat_u, start = 10))
  expect_error(life_expect(mat_u, start = "stage name"))
  expect_error(life_expect(mat_u_na))
  expect_error(life_expect(mat_u_named, start = "invalid stage"))
  expect_error(life_expect(mat_u_named_mismatch, start = "sm"))
  expect_error(life_expect(mat_u_named_partial, start = "sm"))
})
