
test_that("standard_stages works correctly", {
  
  ## mpm without inter-reproductive stage
  repro_stages <- apply(mat_f, 2, function(x) any(x > 0))
  matrix_stages <- c('active', 'active', 'active', 'active')
  r <- mpm_rearrange(mat_u, mat_f, mat_c, repro_stages, matrix_stages)

  x <- standard_stages(r$matF, r$repro_stages, r$matrix_stages)
  expect_type(x, "list")
  expect_length(x, 4)
  
  # make sure all stages in original mpm represented in output
  stages <- as.numeric(unlist(x))
  stages <- stages[!is.na(stages)]
  expect_true(all(1:nrow(r$matU) %in% stages))

  
  ## mpm with inter-reproductive stage
  repro_stages <- apply(mat_f_inter, 2, function(x) any(x > 0))
  matrix_stages <- c('prop', 'active', 'active', 'active', 'active')
  r <- mpm_rearrange(mat_u_inter, mat_f_inter, repro_stages = repro_stages,
                     matrix_stages = matrix_stages)
  
  x <- standard_stages(r$matF, r$repro_stages, r$matrix_stages)
  expect_type(x, "list")
  expect_length(x, 4)
  
  # make sure all stages in original mpm represented in output
  stages <- as.numeric(unlist(x))
  stages <- stages[!is.na(stages)]
  expect_true(all(1:nrow(r$matU) %in% stages))
  
  
  ## mpm with all stages reproductive
  repro_stages <- apply(mat_f_allrep, 2, function(x) any(x > 0))
  matrix_stages <- c('active', 'active')
  r <- mpm_rearrange(mat_u_allrep, mat_f_allrep, repro_stages = repro_stages,
                     matrix_stages = matrix_stages)
  
  x <- standard_stages(r$matF, r$repro_stages, r$matrix_stages)
  expect_type(x, "list")
  expect_length(x, 4)
  expect_true(is.na(x$propStages))
  expect_true(is.na(x$preRepStages))
  expect_true(is.na(x$postRepStages))

  # make sure all stages in original mpm represented in output
  expect_true(all(1:nrow(r$matU) %in% stages))
})


test_that("standard_stages warns and fails gracefully", {
  
  # arguments of different dimension
  repro_stages <- apply(mat_f, 2, function(x) any(x > 0))
  matrix_stages <- c('prop', 'active', 'active')
  expect_error(standard_stages(mat_f, repro_stages, matrix_stages))
  
  # matF contains NA
  repro_stages <- apply(mat_f_na, 2, function(x) any(x > 0))
  matrix_stages <- c('prop', 'active', 'active', 'active')
  expect_error(standard_stages(mat_f_na, repro_stages, matrix_stages))
  
  # no reproductive stages
  repro_stages <- apply(mat_f_zero, 2, function(x) any(x > 0))
  matrix_stages <- c('active', 'active', 'active', 'active')
  expect_error(standard_stages(mat_f_zero, repro_stages, matrix_stages))
})
