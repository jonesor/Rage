context("standardized_stages")

test_that("standardized_stages works correctly", {
  
  ## mpm without inter-reproductive stage
  reproStages <- apply(mat_f, 2, function(x) any(x > 0))
  matrixStages <- c('active', 'active', 'active', 'active')
  r <- mpm_rearrange(mat_u, mat_f, mat_c, reproStages, matrixStages)

  x <- standardized_stages(r$matF, r$reproStages, r$matrixStages)
  expect_is(x, "list")
  expect_length(x, 4)
  
  # make sure all stages in original mpm represented in output
  stages <- as.numeric(unlist(x))
  stages <- stages[!is.na(stages)]
  expect_true(all(stages %in% 1:nrow(r$matU)))

  
  ## mpm with inter-reproductive stage
  reproStages <- apply(mat_f_inter, 2, function(x) any(x > 0))
  matrixStages <- c('prop', 'active', 'active', 'active', 'active')
  r <- mpm_rearrange(mat_u_inter, mat_f_inter, reproStages = reproStages,
                     matrixStages = matrixStages)
  
  x <- standardized_stages(r$matF, r$reproStages, r$matrixStages)
  expect_is(x, "list")
  expect_length(x, 4)
  
  # make sure all stages in original mpm represented in output
  stages <- as.numeric(unlist(x))
  stages <- stages[!is.na(stages)]
  expect_true(all(stages %in% 1:nrow(r$matU)))
  
  
  ## mpm with all stages reproductive
  reproStages <- apply(mat_f_allrep, 2, function(x) any(x > 0))
  matrixStages <- c('active', 'active')
  r <- mpm_rearrange(mat_u_allrep, mat_f_allrep, reproStages = reproStages,
                     matrixStages = matrixStages)
  
  x <- standardized_stages(r$matF, r$reproStages, r$matrixStages)
  expect_is(x, "list")
  expect_length(x, 4)
  expect_true(is.na(x$propStages))
  expect_true(is.na(x$preRepStages))
  expect_true(is.na(x$postRepStages))
})


test_that("standardized_stages warns and fails gracefully", {
  
  # arguments of different dimension
  reproStages <- apply(mat_f, 2, function(x) any(x > 0))
  matrixStages <- c('prop', 'active', 'active')
  expect_error(standardized_stages(mat_f, mat_c, reproStages, matrixStages))
  
  # matF contains NA
  reproStages <- apply(mat_f, 2, function(x) any(x > 0))
  matrixStages <- c('prop', 'active', 'active', 'active')
  expect_error(standardized_stages(mat_f_na, reproStages, matrixStages))
  
  # no reproductive stages
  reproStages <- apply(mat_f_zero, 2, function(x) any(x > 0))
  matrixStages <- c('active', 'active', 'active', 'active')
  expect_error(standardized_stages(mat_f_zero, reproStages, matrixStages))
})
