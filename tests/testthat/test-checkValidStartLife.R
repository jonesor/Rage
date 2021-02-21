
test_that("checkValidStartLife works correctly", {
  
  M <- matrix(1:4, nrow = 2)
  
  startLife <- 1
  expect_silent(checkValidStartLife(startLife, M))
  
  # length(startLife) > 1
  startLife <- c(3, 4)
  expect_error(checkValidStartLife(startLife, M))
  
  # outside 1:nrow(M)
  startLife <- 5
  expect_error(checkValidStartLife(startLife, M))
  
  # start_vec = TRUE
  startLife <- c(0.4, 0.6)
  expect_silent(checkValidStartLife(startLife, M, start_vec = TRUE))
  
  startLife <- c(0, 0.4, 0.6)
  expect_error(checkValidStartLife(startLife, M, start_vec = TRUE))
})
