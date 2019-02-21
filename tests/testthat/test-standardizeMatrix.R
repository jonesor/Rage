context("standardizeMatrix")

test_that("standardizeMatrix works correctly", {

  # mpm with inter-reproductive stage
  reproStages <- apply(mat_f_inter, 2, function(x) any(x > 0))
  matrixStages <- c('prop', 'active', 'active', 'active', 'active')

  x <- standardizeMatrix(mat_u_inter, mat_f_inter, reproStages = reproStages,
                         matrixStages = matrixStages)
  
  expect_is(x, "list")
  expect_length(x, 4)
  expect_equal(ncol(x$matA), 4)
})
