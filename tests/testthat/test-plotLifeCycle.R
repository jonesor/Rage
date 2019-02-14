context("plotLifeCycle")

test_that("plotLifeCycle works correctly", {

  mat_a <- mat_u + mat_f
  stages <- letters[1:nrow(mat_a)]
  
  p1 <- plotLifeCycle(mat_a)
  p2 <- plotLifeCycle(mat_a, stages = stages)
  p3 <- plotLifeCycle(mat_a, stages = stages, title = "my_title")
  
  expect_is(p1, "grViz")
  expect_is(p2, "grViz")
  expect_is(p3, "grViz")
})

