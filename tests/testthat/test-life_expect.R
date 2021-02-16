
test_that("life_expect works correctly", {
  
  l0 <- life_expect(mat_u)
  l0_singular <- life_expect(mat_u_singular)

  expect_equal(ncol(l0), 2L)
  expect_true(all(l0 > 0))
  expect_true(all(unlist(is.na(l0_singular))))
})

test_that("life_expect warns and fails gracefully", {
  
  expect_warning(life_expect(mat_u_zero))
  expect_warning(life_expect(mat_u_survissue))
  expect_error(life_expect(mat_u, start = 10))
  expect_error(life_expect(mat_u_na))
})
