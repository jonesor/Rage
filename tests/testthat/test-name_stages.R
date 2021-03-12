
test_that("name_stages works correctly", {
  s <- 10 # matrix dimension
  
  # generate lists of matrices
  x1 <- replicate(2, matrix(runif(s^2), s, s), simplify = FALSE)
  
  y1 <- name_stages(x1)
  expect_true(inherits(y1, "list"))
  expect_true(inherits(y1[[1]], "matrix"))
  expect_true(nrow(y1[[1]]) == s)
  
  y2 <- name_stages(x1[[1]])
  expect_true(inherits(y2, "matrix"))
  expect_true(nrow(y2) == s)
  
  y3 <- name_stages(x1[[1]], prefix = "", left_pad = FALSE)
  expect_true(identical(rownames(y3), as.character(1:s)))
  
  y4 <- name_stages(x1[[1]], names = LETTERS[1:s], prefix = NULL)
  expect_true(identical(rownames(y4), LETTERS[1:s]))
})

test_that("name_stages warns and fails gracefully", {
  s <- 10 # matrix dimension
  sp1 <- s + 1
  
  # generate lists of matrices
  x1 <- x2 <- replicate(2, matrix(runif(s^2), s, s), simplify = FALSE)
  x2[[3]] <- matrix(runif(sp1^2), sp1, sp1)
  
  # passing wrong object types
  expect_error(name_stages("a"))
  expect_error(name_stages(matrix(1:3)))
  expect_error(name_stages(list("a", "b", "c")))
  
  # passing list of matrices with different dimensions
  expect_error(name_stages(x2))
  
  # overwriting existing stage names
  expect_warning(name_stages(mat_u_named))
  
  # supplying both prefix and names arguments
  expect_warning(name_stages(x1, names = LETTERS[1:s]))
  
  # incorrect number of stage names
  expect_error(name_stages(x1, names = letters[1:(s-1)], prefix = NULL))
})
