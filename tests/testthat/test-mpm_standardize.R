test_that("mpm_standardize works correctly", {
  # mpm with inter-reproductive stage
  repro_stages <- apply(mat_f_inter, 2, function(x) any(x > 0))
  matrix_stages <- c("prop", "active", "active", "active", "active")

  x1 <- mpm_standardize(mat_u_inter, mat_f_inter,
    repro_stages = repro_stages,
    matrix_stages = matrix_stages
  )

  expect_type(x1, "list")
  expect_length(x1, 4)
  expect_equal(ncol(x1$matA), 4)

  x2 <- mpm_standardise(mat_u_inter, mat_f_inter,
    repro_stages = repro_stages,
    matrix_stages = matrix_stages
  )
  expect_identical(x1, x2) # check aliasing to GB spelling
})
