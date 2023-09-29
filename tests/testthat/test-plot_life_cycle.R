test_that("plot_life_cycle works correctly", {
  mat_a <- mat_u + mat_f
  stages <- letters[seq_len(nrow(mat_a))]

  p1 <- plot_life_cycle(mat_a)
  p2 <- plot_life_cycle(mat_a, stages = stages)
  p3 <- plot_life_cycle(mat_a, stages = stages, title = "my_title")
  p4 <- plot_life_cycle(mat_a, node_order = 1:4)

  expect_s3_class(p1, "grViz")
  expect_s3_class(p2, "grViz")
  expect_s3_class(p3, "grViz")
  expect_s3_class(p4, "grViz")
})


test_that("plot_life_cycle warns and fails gracefully", {
  mat_a <- mat_u + mat_f
  colnames(mat_a) <- paste0("A", seq_len(nrow(mat_a)))
  rownames(mat_a) <- paste0("B", seq_len(nrow(mat_a)))
  # make rownames different from colnames
  stages <- dimnames(mat_a)[[1]]

  expect_message(plot_life_cycle(mat_a))
})

# Test failure with incorrect dimensions.
mat_a <- mat_u + mat_f

# wrong number of stages
expect_error(plot_life_cycle(mat_a, stages = c("A", "B", "C")))

# wrong number of entries in node_order
expect_error(plot_life_cycle(mat_a, node_order = 1:3))
