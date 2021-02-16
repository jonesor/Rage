
test_that("mpm_rearrange works correctly", {

  # mpm with inter-reproductive stage
  reproStages <- apply(mat_f_inter, 2, function(x) any(x > 0))
  matrixStages <- c('prop', 'active', 'active', 'active', 'active')

  x <- mpm_rearrange(mat_u_inter, mat_f_inter, reproStages = reproStages,
                     matrixStages = matrixStages)
  
  expect_type(x, "list")
  expect_length(x, 6)
  
  expect_equal(sum(mat_u_inter), sum(x$matU))
  expect_equal(sum(mat_f_inter), sum(x$matF))
  expect_equal(sort(matrixStages), sort(x$matrixStages))
  expect_true(x$nonRepInterRep %in% seq_along(reproStages))
  
  
  # mpm without inter-reproductive stage
  reproStages <- apply(mat_f, 2, function(x) any(x > 0))
  matrixStages <- c('prop', 'active', 'active', 'active')
  
  x <- mpm_rearrange(mat_u, mat_f, matC = NULL, reproStages, matrixStages)
  
  expect_equal(mat_u, x$matU)
  expect_equal(mat_f, x$matF)
  expect_equal(matrixStages, x$matrixStages)
  expect_true(is.na(x$nonRepInterRep))
  
  # mpm with no reproduction
  reproStages <- apply(mat_f_zero, 2, function(x) any(x > 0))
  matrixStages <- c('prop', 'active', 'active', 'active')
  
  x <- suppressWarnings(
    mpm_rearrange(mat_u, mat_f_zero, matC = NULL, reproStages, matrixStages)
  )
  
  expect_equal(mat_f_zero, x$matF)
  expect_equal(mat_u, x$matU)
})


test_that("mpm_rearrange warns and fails gracefully", {
  
  reproStages <- apply(mat_f, 2, function(x) any(x > 0))
  matrixStages <- c('prop', 'active', 'active')
  
  expect_error(mpm_rearrange(mat_u, mat_f, reproStages = reproStages,
                             matrixStages = matrixStages))
})
