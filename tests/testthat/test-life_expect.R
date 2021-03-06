
test_that("life_expect_mean works correctly", {
  
  l0 <- life_expect_mean(mat_u)
  l0_singular <- life_expect_mean(mat_u_singular)
  l1 <- life_expect_mean(mat_u_named, start = "sm")

  expect_true(l0 > 0)
  expect_true(is.na(l0_singular))
  expect_equal(l0, l1)
})

test_that("life_expect_mean warns and fails gracefully", {
  
  expect_warning(life_expect_mean(mat_u_zero))
  expect_warning(life_expect_mean(mat_u_survissue))
  expect_error(life_expect_mean(mat_u, start = 10))
  expect_error(life_expect_mean(mat_u, start = "stage name"))
  expect_error(life_expect_mean(mat_u_na))
  expect_error(life_expect_mean(mat_u_named, start = "invalid stage"))
  expect_error(life_expect_mean(mat_u_named_mismatch, start = "sm"))
  expect_error(life_expect_mean(mat_u_named_partial, start = "sm"))
})


test_that("life_expect_var works correctly", {
  
  l0 <- life_expect_var(mat_u)
  l0_singular <- life_expect_var(mat_u_singular)
  l1 <- life_expect_var(mat_u_named, start = "sm")
  
  expect_true(l0 > 0)
  expect_true(is.na(l0_singular))
  expect_equal(l0, l1)
})

test_that("life_expect_var warns and fails gracefully", {
  
  expect_warning(life_expect_var(mat_u_zero))
  expect_warning(life_expect_var(mat_u_survissue))
  expect_error(life_expect_var(mat_u, start = 10))
  expect_error(life_expect_var(mat_u, start = "stage name"))
  expect_error(life_expect_var(mat_u_na))
  expect_error(life_expect_var(mat_u_named, start = "invalid stage"))
  expect_error(life_expect_var(mat_u_named_mismatch, start = "sm"))
  expect_error(life_expect_var(mat_u_named_partial, start = "sm"))
})
