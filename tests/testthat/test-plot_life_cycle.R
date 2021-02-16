
test_that("plot_life_cycle works correctly", {

  mat_a <- mat_u + mat_f
  stages <- letters[1:nrow(mat_a)]
  
  p1 <- plot_life_cycle(mat_a)
  p2 <- plot_life_cycle(mat_a, stages = stages)
  p3 <- plot_life_cycle(mat_a, stages = stages, title = "my_title")
  
  expect_s3_class(p1, "grViz")
  expect_s3_class(p2, "grViz")
  expect_s3_class(p3, "grViz")
})

