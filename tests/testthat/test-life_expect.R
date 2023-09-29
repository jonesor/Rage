test_that("life_expect_mean works correctly", {
  l0 <- life_expect_mean(mat_u)
  l0_singular <- life_expect_mean(mat_u_singular)
  l1 <- life_expect_mean(mat_u_named, start = "sm")

  l0_1 <- life_expect_mean(mat_u, start = NULL)
  l0_2 <- life_expect_mean(mat_u, start = NULL, mixdist = c(1, 2, 3, 4))
  l0_3 <- life_expect_mean(mat_u, start = 3)

  expect_gt(l0, 0)
  expect_true(is.na(l0_singular))
  expect_identical(l0, l1)
  expect_vector(l0_1)
  expect_vector(l0_2)
  expect_vector(l0_3)

  l0_s2 <- life_expect_mean(mat_u, start = 2)
  l0_s2_mix <- life_expect_mean(mat_u, start = NULL, mixdist = c(0, 1, 0, 0))
  expect_identical(l0_s2, l0_s2_mix)
})

test_that("life_expect_mean warns and fails gracefully", {
  expect_warning(life_expect_mean(mat_u_zero))
  expect_warning(life_expect_mean(mat_u_survissue))
  expect_error(life_expect_mean(mat_u, start = 10))
  expect_error(life_expect_mean(mat_u, start = "stage name"))
  expect_error(life_expect_mean(mat_u_na))
  expect_error(life_expect_mean(mat_u_named, start = "invalid stage"))
  expect_error(life_expect_mean(mat_u, start = 2, mixdist = c(1, 2, 3, 4)))
  expect_error(life_expect_mean(mat_u, start = 10))
  expect_error(life_expect_mean(mat_u, start = 0))
})


test_that("life_expect_var works correctly", {
  l0 <- life_expect_var(mat_u)
  l0_singular <- life_expect_var(mat_u_singular)
  l1 <- life_expect_var(mat_u_named, start = "sm")

  expect_gt(l0, 0)
  expect_true(is.na(l0_singular))
  expect_identical(l0, l1)

  l0_1 <- life_expect_var(mat_u, start = NULL)
  l0_2 <- life_expect_mean(mat_u, start = NULL, mixdist = c(1, 2, 3, 4))
  l0_3 <- life_expect_mean(mat_u, start = 3)
  expect_vector(l0_1)
  expect_vector(l0_2)
  expect_vector(l0_3)
})

test_that("life_expect_var warns and fails gracefully", {
  expect_warning(life_expect_var(mat_u_zero))
  expect_warning(life_expect_var(mat_u_survissue))
  expect_error(life_expect_var(mat_u, start = 10))
  expect_error(life_expect_var(mat_u, start = "stage name"))
  expect_error(life_expect_var(mat_u_na))
  expect_error(life_expect_var(mat_u_named, start = "invalid stage"))
})
