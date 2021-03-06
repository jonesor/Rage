
test_that("nameStages works correctly", {
  s <- 10 # matrix dimension
  
  # generate lists of matrices
  x1 <- replicate(2, matrix(runif(s^2), s, s), simplify = FALSE)
  
  y1 <- nameStages(x1)
  expect_true(inherits(y1, "list"))
  expect_true(inherits(y1[[1]], "matrix"))
  expect_true(nrow(y1[[1]]) == s)
  
  y2 <- nameStages(x1[[1]])
  expect_true(inherits(y2, "matrix"))
  expect_true(nrow(y2) == s)
  
  y3 <- nameStages(x1[[1]], prefix = "", left_pad = FALSE)
  expect_true(identical(rownames(y3), as.character(1:s)))
  
  y4 <- nameStages(x1[[1]], names = LETTERS[1:s], prefix = NULL)
  expect_true(identical(rownames(y4), LETTERS[1:s]))
})

test_that("nameStages warns and fails gracefully", {
  s <- 10 # matrix dimension
  sp1 <- s + 1
  
  # generate lists of matrices
  x1 <- x2 <- replicate(2, matrix(runif(s^2), s, s), simplify = FALSE)
  x2[[3]] <- matrix(runif(sp1^2), sp1, sp1)
  
  # passing wrong object types
  expect_error(nameStages("a"))
  expect_error(nameStages(matrix(1:3)))
  expect_error(nameStages(list("a", "b", "c")))
  
  # passing list of matrices with different dimensions
  expect_error(nameStages(x2))
  
  # overwriting existing stage names
  expect_warning(nameStages(x3))
  
  # supplying both prefix and names arguments
  expect_warning(nameStages(x1, names = LETTERS[1:s]))
  
  # incorrect number of stage names
  expect_error(nameStages(x1, names = letters[1:(s-1)], prefix = NULL))
})
