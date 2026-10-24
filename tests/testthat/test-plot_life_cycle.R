# ---- graphviz engine ---------------------------------------------------------

test_that("plot_life_cycle graphviz engine works correctly", {
  skip_if_not_installed("DiagrammeR")
  mat_a <- mat_u + mat_f
  stages <- letters[seq_len(nrow(mat_a))]
  mat_a_na <- mat_a
  mat_a_na[1, 3] <- NA

  p1 <- plot_life_cycle(mat_a, engine = "graphviz")
  p2 <- plot_life_cycle(mat_a, stages = stages, engine = "graphviz")
  p3 <- plot_life_cycle(mat_a, stages = stages, title = "my_title",
                        engine = "graphviz")
  p4 <- plot_life_cycle(mat_a, node_order = 1:4, engine = "graphviz")
  p5 <- plot_life_cycle(mat_a_na, engine = "graphviz")
  p6 <- plot_life_cycle(matU = mat_u, matF = mat_f, engine = "graphviz")

  expect_s3_class(p1, "grViz")
  expect_s3_class(p2, "grViz")
  expect_s3_class(p3, "grViz")
  expect_s3_class(p4, "grViz")
  expect_s3_class(p5, "grViz")
  expect_s3_class(p6, "grViz")
  expect_identical(p1$x$diagram, p6$x$diagram)
  expect_match(p1$x$diagram, "graph\\[overlap=false, margin=0, pad=0\\]")
  expect_false(grepl("x \\[style=invis\\]", p1$x$diagram))
  expect_false(grepl("labelloc=\"t\"", p1$x$diagram, fixed = TRUE))
  expect_true(grepl("labelloc=\"t\"", p3$x$diagram, fixed = TRUE))
  expect_true(grepl("\n  label=\"my_title\"", p3$x$diagram, fixed = TRUE))
})

test_that("plot_life_cycle graphviz engine warns and fails gracefully", {
  skip_if_not_installed("DiagrammeR")
  mat_a <- mat_u + mat_f
  colnames(mat_a) <- paste0("A", seq_len(nrow(mat_a)))
  rownames(mat_a) <- paste0("B", seq_len(nrow(mat_a)))

  expect_message(plot_life_cycle(mat_a, engine = "graphviz"))
  expect_error(plot_life_cycle(mat_u + mat_f, stages = c("A", "B", "C"),
                               engine = "graphviz"))
  expect_error(plot_life_cycle(mat_u + mat_f, node_order = 1:3,
                               engine = "graphviz"))
})

# ---- ggplot engine (default) ---------------------------------------------------

test_that("plot_life_cycle ggplot engine returns a sized ggplot", {
  mat_a <- mat_u + mat_f
  p <- plot_life_cycle(mat_a)

  expect_s3_class(p, "ggplot")
  sz <- attr(p, "lc_size")
  expect_named(sz, c("width", "height"))
  expect_true(all(sz > 0))

  # builds without error or warning
  expect_silent(ggplot2::ggplot_build(p))

  # titles, stage names and node order are accepted
  p2 <- plot_life_cycle(mat_a, stages = letters[1:4], title = "my_title",
                        node_order = 4:1)
  expect_s3_class(p2, "ggplot")
  expect_identical(p2$labels$title, "my_title")
  expect_gt(attr(p2, "lc_size")[["height"]], 0)
})

test_that("plot_life_cycle ggplot engine distinguishes transition types", {
  p_a <- plot_life_cycle(mat_u + mat_f)
  p_uf <- plot_life_cycle(matU = mat_u, matF = mat_f)

  types_a <- unique(p_a$layers[[1]]$data$type)
  types_uf <- unique(p_uf$layers[[1]]$data$type)
  expect_identical(types_a, "U")
  expect_setequal(types_uf, c("U", "F"))

  # one arrow per non-zero entry in each matrix
  n_paths <- length(unique(p_uf$layers[[1]]$data$id))
  expect_equal(n_paths, sum(mat_u != 0) + sum(mat_f != 0))

  # legend only when there is more than one transition type
  expect_identical(p_a$theme$legend.position, "none")
  expect_identical(p_uf$theme$legend.position, "bottom")

  mat_c <- mat_f * 0
  mat_c[2, 4] <- 0.3
  p_ufc <- plot_life_cycle(matU = mat_u, matF = mat_f, matC = mat_c)
  expect_setequal(unique(p_ufc$layers[[1]]$data$type), c("U", "F", "C"))
})

test_that("plot_life_cycle ggplot engine labels values and NAs", {
  mat_a <- mat_u + mat_f
  mat_a[1, 3] <- NA
  p <- plot_life_cycle(mat_a, digits = 1)

  labs <- p$layers[[length(p$layers)]]$data$label
  expect_true("NA" %in% labs)
  expect_false(any(grepl("\\.\\d{2,}", labs)))

  paths <- p$layers[[1]]$data
  expect_true(any(paths$na))

  p_novals <- plot_life_cycle(mat_a, show_values = FALSE)
  expect_length(p_novals$layers, length(p$layers) - 1)

  p_scaled <- plot_life_cycle(mat_a, scale_width = TRUE)
  expect_gt(length(unique(p_scaled$layers[[1]]$data$lw)), 1)
})

test_that("plot_life_cycle ggplot engine warns and fails gracefully", {
  mat_a <- mat_u + mat_f
  bad <- mat_a
  colnames(bad) <- paste0("A", seq_len(nrow(bad)))
  rownames(bad) <- paste0("B", seq_len(nrow(bad)))

  expect_message(plot_life_cycle(bad))
  expect_error(plot_life_cycle(mat_a, stages = c("A", "B", "C")))
  expect_error(plot_life_cycle(mat_a, node_order = 1:3))
  expect_error(plot_life_cycle(mat_a, node_order = c(1, 1, 2, 3)))
  expect_error(plot_life_cycle(mat_a, matU = mat_u))
  expect_error(plot_life_cycle(mat_a, matF = mat_f))
  expect_error(plot_life_cycle())
  expect_error(plot_life_cycle(matU = mat_u, matF = mat_f[1:3, 1:3]))
  expect_error(plot_life_cycle(mat_a * 0))
})

test_that("save_life_cycle writes a file at the natural size", {
  p <- plot_life_cycle(matU = mat_u, matF = mat_f)
  f <- tempfile(fileext = ".pdf")
  on.exit(unlink(f))
  expect_identical(save_life_cycle(p, f), f)
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)

  expect_error(save_life_cycle(ggplot2::ggplot(), f), "no stored size")
})

test_that("plot_life_cycle ggplot engine output is visually stable", {
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  mat_a <- mat_u + mat_f
  mat_a_na <- mat_a
  mat_a_na[1, 3] <- NA
  vdiffr::expect_doppelganger(
    "life cycle matA",
    plot_life_cycle(mat_a, base_family = "sans")
  )
  vdiffr::expect_doppelganger(
    "life cycle matU matF named",
    plot_life_cycle(matU = mat_u, matF = mat_f, stages = c(
      "Seed", "Seedling", "Juvenile", "Adult"
    ), title = "Example", base_family = "sans")
  )
})
