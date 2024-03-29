test_that("life_elas works correctly", {
  lx <- 0.8^(0:20)

  x <- life_elas(lx)
  x_trap <- life_elas(lx, trapeze = TRUE)

  expect_length(x, 1L)
  expect_length(x_trap, 1L)
  expect_gte(x, 0)
  expect_gte(x_trap, 0)
})

test_that("life_elas warns and fails gracefully", {
  lx1 <- c(1.1, 0.6, 0.5, 0.4) # lx is greater than 1
  expect_error(life_elas(lx1))

  lx2 <- c(1.0, 0.6, 0.61, 0.5) # lx not monotonically declining
  expect_error(life_elas(lx2))
})
