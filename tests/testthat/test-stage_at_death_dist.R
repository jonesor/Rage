test_that("stage_at_death_dist works correctly", {
  x_default <- stage_at_death_dist(mat_u)
  x_start <- stage_at_death_dist(mat_u, start = c(1, 0, 0, 0))
  x_named <- stage_at_death_dist(mat_u_named)

  expect_type(x_default, "double")
  expect_length(x_default, ncol(mat_u))
  expect_equal(sum(x_default), 1)

  expect_type(x_start, "double")
  expect_length(x_start, ncol(mat_u))
  expect_equal(sum(x_start), 1)
  expect_false(isTRUE(all.equal(x_default, x_start)))

  expect_identical(names(x_named), rownames(mat_u_named))
})

test_that("stage_at_death_dist warns and fails gracefully", {
  expect_error(stage_at_death_dist(mat_u_notsq))
  expect_error(stage_at_death_dist(mat_u, start = "stage_1"))
  expect_error(stage_at_death_dist(mat_u, start = c(1, 0, 0)))
  expect_error(stage_at_death_dist(mat_u, start = c(1, -1, 0, 0)))
  expect_error(stage_at_death_dist(mat_u, start = c(0, 0, 0, 0)))
  expect_error(stage_at_death_dist(mat_u_survissue))
})
