#' Life cycle diagram drawn with ggplot2 (internal engine for plot_life_cycle)
#'
#' @param mats Named list of matrices ("U", and optionally "F", "C"), already
#'   validated as square and of equal dimension.
#' @noRd
.plot_life_cycle_gg <- function(mats, stages, node_order = NULL, title = NULL,
                                digits = 2, show_values = TRUE,
                                scale_width = FALSE,
                                colours = c(U = "#3D4F5C", F = "#D55E00",
                                            C = "#009E73"),
                                type_labels = c(U = "Survival / growth",
                                                F = "Sexual reproduction",
                                                C = "Clonal reproduction"),
                                node_fill = "white", node_line = "#3D4F5C",
                                base_size = 11, base_family = "") {
  n <- nrow(mats[[1]])
  stages <- as.character(stages)
  if (length(mats) == 1) type_labels["U"] <- "Transition"

  if (is.null(node_order)) node_order <- seq_len(n)
  if (length(node_order) != n || !setequal(node_order, seq_len(n))) {
    stop("node_order must be a permutation of 1:", n, call. = FALSE)
  }
  # As in the graphviz engine, stages are displayed in the order
  # stages[order(node_order)], i.e. node_order[i] is the position of stage i
  pos <- as.integer(node_order)

  # ---- Geometry -------------------------------------------------------------
  wrapped <- vapply(stages, function(s) paste(strwrap(s, 10), collapse = "\n"),
                    character(1))
  # Size nodes from the physical size of their labels. One data unit is
  # drawn as `unit_in` inches, so text (in pt) and geometry stay in
  # proportion when the plot is saved at its natural size (see
  # save_life_cycle()).
  unit_in <- base_size / 11
  # measure on a temporary null device, then restore the user's device
  old_dev <- grDevices::dev.cur()
  grDevices::pdf(NULL)
  txt_w <- vapply(wrapped, function(s) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(
      s, gp = grid::gpar(fontsize = base_size, fontface = "bold",
                         fontfamily = base_family, lineheight = 0.85)
    )), "inches", valueOnly = TRUE)
  }, numeric(1))
  txt_h <- vapply(wrapped, function(s) {
    grid::convertHeight(grid::grobHeight(grid::textGrob(
      s, gp = grid::gpar(fontsize = base_size, fontface = "bold",
                         fontfamily = base_family, lineheight = 0.85)
    )), "inches", valueOnly = TRUE)
  }, numeric(1))
  grDevices::dev.off()
  if (old_dev > 1) grDevices::dev.set(old_dev)
  R <- max(0.28, (sqrt(txt_w^2 + txt_h^2) / 2) / unit_in + 0.09)
  gap <- max(1.15, 2 * R + 0.75)            # centre-to-centre distance
  loop_r <- 0.20                             # self-loop radius
  node_x <- (pos - 1) * gap

  circle_pts <- function(cx, cy, r, from = 0, to = 2 * pi, k = 80) {
    a <- seq(from, to, length.out = k)
    data.frame(x = cx + r * cos(a), y = cy + r * sin(a))
  }

  # ---- Build edge list ------------------------------------------------------
  edges <- do.call(rbind, lapply(names(mats), function(nm) {
    m <- mats[[nm]]
    idx <- which(is.na(m) | m != 0, arr.ind = TRUE)
    if (nrow(idx) == 0) return(NULL)
    data.frame(from = idx[, "col"], to = idx[, "row"],
               value = m[idx], type = nm, stringsAsFactors = FALSE)
  }))
  if (is.null(edges)) stop("All matrix entries are zero; nothing to draw", call. = FALSE)

  edges$label <- vapply(seq_len(nrow(edges)), function(i) {
    format(round(edges$value[i], digits), drop0trailing = TRUE,
           scientific = FALSE, nsmall = 0)
  }, character(1))
  edges$label[is.na(edges$value)] <- "NA"

  # Width: optional scaling within each type
  lw0 <- base_size / 22
  edges$lw <- lw0
  if (scale_width) {
    for (tp in unique(edges$type)) {
      i <- edges$type == tp & !is.na(edges$value)
      v <- edges$value[i]
      edges$lw[i] <- lw0 * (0.5 + 2 * v / max(v))
    }
  }

  # ---- Classify edges and assign "ports" on node boundaries -------------
  p_from <- pos[edges$from]; p_to <- pos[edges$to]
  edges$span <- abs(p_to - p_from)
  edges$kind <- ifelse(p_from == p_to, "loop",
                ifelse(edges$span == 1 & p_to > p_from & edges$type == "U",
                       "straight", "arc"))
  edges$up <- ifelse(edges$kind == "loop",
                     ifelse(edges$type == "U", 1, -1),
                     ifelse(p_to > p_from & edges$type == "U", 1, -1))
  # stacking index for arcs joining the same pair of positions on one side
  key <- paste(pmin(p_from, p_to), pmax(p_from, p_to), edges$up)
  edges$k <- stats::ave(seq_along(key), key, FUN = seq_along) - 1
  edges$k[edges$kind != "arc"] <- 0
  # Depth level per side: an arc sits one level below everything it encloses,
  # and is pushed further down if its label would collide with an
  # overlapping arc already at that level
  edges$level <- 0
  for (sd in c(-1, 1)) {
    ii <- which(edges$kind == "arc" & edges$up == sd)
    ii <- ii[order(edges$span[ii])]
    for (i in ii) {
      a <- min(p_from[i], p_to[i]); b <- max(p_from[i], p_to[i])
      done <- ii[edges$level[ii] > 0 & ii != i]
      lo <- pmin(p_from[done], p_to[done]); hi <- pmax(p_from[done], p_to[done])
      inside <- done[lo >= a & hi <= b & !(lo == a & hi == b)]
      lev <- if (length(inside)) max(edges$level[inside]) + 1 else 1
      same <- done[lo == a & hi == b]
      if (length(same)) { edges$level[i] <- edges$level[same[1]]; next }
      repeat {
        clash <- done[edges$level[done] == lev &
                        abs((lo + hi) / 2 - (a + b) / 2) < 1.5 &
                        lo < b & hi > a]
        if (!length(clash)) break
        lev <- lev + 1
      }
      edges$level[i] <- lev
    }
  }

  # Each arc has two endpoints. On each node side, endpoints heading right
  # are spread over the right quadrant and those heading left over the left
  # quadrant, ordered by span (longer, taller arcs nearer the vertical) so
  # arcs never cross at the node and never share an entry point.
  edges$th_from <- NA_real_; edges$th_to <- NA_real_
  arc_i <- which(edges$kind == "arc")
  if (length(arc_i)) {
    ends <- data.frame(
      edge = rep(arc_i, 2),
      end = rep(c("from", "to"), each = length(arc_i)),
      node = c(p_from[arc_i], p_to[arc_i]),
      other = c(p_to[arc_i], p_from[arc_i]),
      up = rep(edges$up[arc_i], 2),
      span = rep(edges$span[arc_i], 2),
      k = rep(edges$k[arc_i], 2)
    )
    ends$right <- ends$other > ends$node
    loop_side <- paste(pos[edges$from], edges$up)[edges$kind == "loop"]
    side_key <- paste(ends$node, ends$up)
    grp <- split(seq_len(nrow(ends)), paste(side_key, ends$right))
    for (g in grp) {
      has_loop <- side_key[g[1]] %in% loop_side
      # if nothing leaves this side in the other direction, use the whole side
      other_dir <- any(side_key == side_key[g[1]] & ends$right != ends$right[g[1]])
      hi <- if (has_loop) 0.27 * pi else if (other_dir) 0.46 * pi else 0.62 * pi
      lo <- 0.12 * pi
      o <- g[order(ends$span[g], ends$k[g])]
      m <- length(o)
      th <- if (m == 1) {
        min(hi, lo + 0.07 * pi * ends$span[o])
      } else {
        lo + (seq_len(m) - 0.5) / m * (hi - lo)
      }
      for (j in seq_len(m)) {
        e <- ends$edge[o[j]]
        if (ends$end[o[j]] == "from") edges$th_from[e] <- th[j]
        else edges$th_to[e] <- th[j]
      }
    }
  }

  # Geometry for each edge: returns path points and a label anchor
  edge_geom <- function(e) {
    k <- e$k
    p1 <- pos[e$from]; p2 <- pos[e$to]
    x1 <- node_x[e$from]; x2 <- node_x[e$to]
    span <- abs(p2 - p1)
    if (p1 == p2) {
      # Self-loop: survival above, reproduction below
      up <- if (e$type == "U") 1 else -1
      cy <- up * (R + loop_r * 0.55)
      ang0 <- if (up > 0) -0.25 * pi else 0.75 * pi
      # trace the loop outside the node, starting and ending on its boundary
      a <- seq(ang0 + 0.55, ang0 + 2 * pi - 0.55, length.out = 60)
      if (up < 0) a <- rev(a)
      pts <- data.frame(x = x1 + loop_r * cos(a), y = cy + loop_r * sin(a))
      pts <- pts[(pts$x - x1)^2 + pts$y^2 >= R^2 * 0.98, ]
      lab <- c(x1, cy + up * (loop_r + 0.13))
    } else if (span == 1 && p2 > p1 && e$type == "U") {
      # Adjacent forward step: straight arrow
      dir <- sign(x2 - x1)
      pts <- data.frame(x = c(x1 + dir * R, x2 - dir * R), y = c(0, 0))
      lab <- c((x1 + x2) / 2, 0.14)
    } else {
      # Arc: forward survival above; everything else below.
      up <- if (p2 > p1 && e$type == "U") 1 else -1
      # stack arcs with the same span/side (e.g. U and F both 2 -> 1)
      h <- R + 0.05 + 0.5 * e$level + 0.22 * k
      if (up > 0) h <- h + loop_r * 2       # clear the stasis loops
      dir <- sign(x2 - x1)
      s <- c(x1 + dir * R * cos(e$th_from), up * R * sin(e$th_from))
      f <- c(x2 - dir * R * cos(e$th_to), up * R * sin(e$th_to))
      # cubic Bezier leaving each node radially, so arcs whose ports are
      # ordered by span nest without crossing; L sets the peak height to h
      # tangent is a blend of the radial and vertical directions: radial
      # enough to look natural, vertical enough not to kink on short spans
      ps <- 0.65 * pi / 2 + 0.35 * e$th_from
      pf <- 0.65 * pi / 2 + 0.35 * e$th_to
      rs <- c(dir * cos(ps), up * sin(ps))
      rf <- c(-dir * cos(pf), up * sin(pf))
      L <- (h - (abs(s[2]) + abs(f[2])) / 2) * 8 /
        (3 * (abs(rs[2]) + abs(rf[2])))
      L <- max(L, 0.1)
      c1 <- s + L * rs; c2 <- f + L * rf
      t <- seq(0, 1, length.out = 99)
      b0 <- (1 - t)^3; b1 <- 3 * (1 - t)^2 * t; b2 <- 3 * (1 - t) * t^2
      b3 <- t^3
      pts <- data.frame(x = b0 * s[1] + b1 * c1[1] + b2 * c2[1] + b3 * f[1],
                        y = b0 * s[2] + b1 * c1[2] + b2 * c2[2] + b3 * f[2])
      lab <- c(pts$x[50], pts$y[50])
    }
    list(pts = pts, lab = lab)
  }

  paths <- vector("list", nrow(edges))
  labs <- vector("list", nrow(edges))
  for (i in seq_len(nrow(edges))) {
    g <- edge_geom(edges[i, ])
    paths[[i]] <- cbind(g$pts, id = i, type = edges$type[i],
                        lw = edges$lw[i], na = is.na(edges$value[i]))
    labs[[i]] <- data.frame(x = g$lab[1], y = g$lab[2],
                            label = edges$label[i], type = edges$type[i])
  }
  paths <- do.call(rbind, paths)
  labs <- do.call(rbind, labs)

  nodes <- do.call(rbind, lapply(seq_len(n), function(i) {
    cbind(circle_pts(node_x[i], 0, R), id = i)
  }))
  node_lab <- data.frame(x = node_x, y = 0, label = wrapped)

  # Panel extents, so the saved size maps one data unit to unit_in inches
  lab_pad <- 0.16
  xr <- range(c(paths$x, labs$x, node_x - R, node_x + R)) + c(-0.15, 0.15)
  yr <- range(c(paths$y, labs$y - lab_pad, labs$y + lab_pad, -R, R)) +
    c(-0.05, 0.05)

  present <- intersect(c("U", "F", "C"), unique(edges$type))
  pt <- base_size / ggplot2::.pt
  arrow_spec <- grid::arrow(type = "closed", angle = 22,
                            length = grid::unit(base_size * 0.22, "mm"))

  # ---- Plot -----------------------------------------------------------------
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = paths,
      ggplot2::aes(x, y, group = id, colour = type, linewidth = lw,
                   linetype = na),
      arrow = arrow_spec, linejoin = "mitre", lineend = "butt"
    ) +
    ggplot2::geom_polygon(
      data = nodes, ggplot2::aes(x, y, group = id),
      fill = node_fill, colour = node_line, linewidth = lw0 * 1.4
    ) +
    ggplot2::geom_text(
      data = node_lab, ggplot2::aes(x, y, label = label),
      size = pt, family = base_family, fontface = "bold",
      colour = node_line, lineheight = 0.85
    ) +
    ggplot2::scale_colour_manual(values = colours[present],
                                 labels = type_labels[present],
                                 breaks = present, name = NULL) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid",
                                              `TRUE` = "22"),
                                   guide = "none") +
    ggplot2::coord_fixed(xlim = xr, ylim = yr, expand = FALSE,
                         clip = "off") +
    ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      legend.position = if (length(present) > 1) "bottom" else "none",
      legend.key.width = grid::unit(base_size * 1.4, "pt"),
      legend.text = ggplot2::element_text(
        margin = ggplot2::margin(l = 3, r = base_size)
      ),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0,
                                         margin = ggplot2::margin(b = 6)),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(linewidth = lw0 * 1.5, arrow = NULL)
    ))

  if (show_values) {
    p <- p + ggplot2::geom_label(
      data = labs, ggplot2::aes(x, y, label = label, colour = type),
      size = pt * 0.82, family = base_family, fill = "white",
      label.size = 0, label.padding = grid::unit(0.12, "lines"),
      show.legend = FALSE
    )
  }
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)

  # Natural size in inches, used by save_life_cycle()
  margin_in <- 2 * 8 / 72
  attr(p, "lc_size") <- c(
    width = diff(xr) * unit_in + margin_in,
    height = diff(yr) * unit_in + margin_in +
      (if (length(present) > 1) base_size * 2.6 / 72 else 0) +
      (if (!is.null(title)) base_size * 2.2 / 72 else 0)
  )
  p
}


#' Save a life cycle diagram at its natural size
#'
#' Saves a diagram made by [plot_life_cycle()] (with `engine = "ggplot"`) using
#' [ggplot2::ggsave()], at the width and height for which text, nodes and
#' arrows are in proportion. The file format follows the extension of
#' `filename`: `.pdf` and `.svg` give vector output suitable for journals and
#' for PowerPoint/Keynote; `.png` gives a raster image.
#'
#' @param plot A diagram returned by [plot_life_cycle()] with
#'   `engine = "ggplot"`.
#' @param filename File name to write.
#' @param scale Multiplier applied to the natural width and height. Text and
#'   lines scale with it.
#' @param ... Further arguments passed to [ggplot2::ggsave()], such as `dpi`
#'   or `device`.
#'
#' @return Invisibly, the file name.
#'
#' @details The natural size is stored as an attribute of the plot. Adding
#'   further ggplot2 layers or themes with `+` drops it; in that case call
#'   [ggplot2::ggsave()] directly with your own `width` and `height`.
#'
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#'
#' @family visualisation
#'
#' @examples
#' \dontrun{
#' data(mpm1)
#' p <- plot_life_cycle(matU = mpm1$matU, matF = mpm1$matF)
#' save_life_cycle(p, "life_cycle.pdf")
#' save_life_cycle(p, "life_cycle.png", dpi = 300)
#' }
#'
#' @export save_life_cycle
save_life_cycle <- function(plot, filename, scale = 1, ...) {
  sz <- attr(plot, "lc_size")
  if (is.null(sz)) {
    stop("'plot' has no stored size. Use a plot made by ",
         "plot_life_cycle(engine = \"ggplot\"), or call ggplot2::ggsave() ",
         "with your own width and height.", call. = FALSE)
  }
  ggplot2::ggsave(filename, plot,
                  width = sz[["width"]] * scale,
                  height = sz[["height"]] * scale,
                  units = "in", bg = "white", ...)
  invisible(filename)
}
