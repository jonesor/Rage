
test_that("mpm_rearrange works correctly", {

  # mpm with inter-reproductive stage
  repro_stages <- apply(mat_f_inter, 2, function(x) any(x > 0))
  matrix_stages <- c("prop", "active", "active", "active", "active")

  x <- mpm_rearrange(mat_u_inter, mat_f_inter,
    repro_stages = repro_stages,
    matrix_stages = matrix_stages
  )

  expect_type(x, "list")
  expect_length(x, 6)

  expect_equal(sum(mat_u_inter), sum(x$matU))
  expect_equal(sum(mat_f_inter), sum(x$matF))
  expect_equal(sort(matrix_stages), sort(x$matrix_stages))
  expect_true(x$nonRepInterRep %in% seq_along(repro_stages))


  # mpm without inter-reproductive stage
  repro_stages <- apply(mat_f, 2, function(x) any(x > 0))
  matrix_stages <- c("prop", "active", "active", "active")

  x <- mpm_rearrange(mat_u, mat_f, matC = NULL, repro_stages, matrix_stages)

  expect_equal(mat_u, x$matU)
  expect_equal(mat_f, x$matF)
  expect_equal(matrix_stages, x$matrix_stages)
  expect_true(is.na(x$nonRepInterRep))

  # mpm with no reproduction
  repro_stages <- apply(mat_f_zero, 2, function(x) any(x > 0))
  matrix_stages <- c("prop", "active", "active", "active")

  x <- suppressWarnings(
    mpm_rearrange(mat_u, mat_f_zero, matC = NULL, repro_stages, matrix_stages)
  )

  expect_equal(mat_f_zero, x$matF)
  expect_equal(mat_u, x$matU)
})


test_that("mpm_rearrange warns and fails gracefully", {
  repro_stages <- apply(mat_f, 2, function(x) any(x > 0))
  matrix_stages <- c("prop", "active", "active")

  expect_error(mpm_rearrange(mat_u, mat_f,
    repro_stages = repro_stages,
    matrix_stages = matrix_stages
  ))
})
