
test_that("perturb_stochastic works correctly", {
  s <- 3 # matrix dimension
  N <- 15 # length of time series

  # generate lists of matrices and population vectors
  X1 <- replicate(N, matrix(runif(s^2), s, s), simplify = FALSE)
  u1 <- replicate(N, runif(s), simplify = FALSE)
  u1 <- lapply(u1, function(x) x / sum(x))

  e1 <- perturb_stochastic(X1, u1)
  expect_type(e1, "list")
  expect_true(inherits(e1[[1]], "matrix"))
  expect_length(e1, 3)
  expect_true(nrow(e1[[1]]) == s)


  ### test failure
  X2 <- X1[1:(N - 1)]
  u2 <- u1
  expect_error(perturb_stochastic(X2, u2))

  X3 <- X1
  u3 <- u1
  X3[[1]] <- X3[[1]][1:(s - 1), 1:(s - 1)]
  expect_error(perturb_stochastic(X3, u3))

  X4 <- X1
  u4 <- u1
  u4[[1]] <- u4[[1]][1:(s - 1)]
  expect_error(perturb_stochastic(X4, u4))

  X5 <- X1
  u5 <- u1
  X5[[1]][1, 1] <- NA
  expect_error(perturb_stochastic(X5, u5))
})
