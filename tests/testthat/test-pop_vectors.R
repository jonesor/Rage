test_that("pop_vectors works correctly", {
  s <- 3 # matrix dimension
  N <- 7 # length of time series

  # generate lists of matrices
  x1 <- replicate(N, matrix(runif(s^2), s, s), simplify = FALSE)

  test_sum_1 <- function(x) abs(sum(x) - 1) < 0.000001

  v1 <- pop_vectors(x1)
  expect_type(v1, "list")
  expect_type(v1[[1]], "double")
  expect_length(v1, N)
  expect_length(v1[[1]], s)
  expect_true(all(vapply(v1, test_sum_1, logical(1))))

  v2 <- pop_vectors(x1, start = "uniform")
  expect_type(v2, "list")
  expect_type(v2[[1]], "double")
  expect_length(v2, N)
  expect_length(v2[[1]], s)
  expect_true(all(vapply(v2, test_sum_1, logical(1))))

  v3 <- pop_vectors(x1, start = "random")
  expect_type(v3, "list")
  expect_type(v3[[1]], "double")
  expect_length(v3, N)
  expect_length(v3[[1]], s)
  expect_true(all(vapply(v3, test_sum_1, logical(1))))
})


test_that("pop_vectors warns and fails gracefully", {
  s <- 3 # matrix dimension
  N <- 7 # length of time series

  x2 <- replicate(N, matrix(runif(s^2), s, s), simplify = FALSE)
  expect_error(pop_vectors(x2, start = "blagjd"))
})
