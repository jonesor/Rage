#' Plot a life cycle diagram from a matrix population model
#'
#' Plots the life cycle diagram illustrated by a matrix population model.
#' Transitions with value \code{0} are assumed to be biologically impossible,
#' and arrows are omitted from the diagram. Transition rates of \code{NA} are
#' shown as arrows labelled \code{NA} (dashed with \code{engine = "ggplot"}).
#'
#' Two drawing engines are available:
#' \describe{
#'   \item{\code{"ggplot"} (default)}{Draws the diagram with \pkg{ggplot2}
#'   using a fixed layout: stages sit in a row, transitions to the next stage
#'   are straight arrows, longer forward transitions arc above, backward
#'   transitions (reproduction, retrogression) arc below, and stasis is a
#'   loop above each stage. If \code{matU} is given together with
#'   \code{matF} and/or \code{matC}, survival/growth, sexual reproduction and
#'   clonal reproduction are distinguished by colour. The result is a ggplot
#'   object, so it can be modified, combined with other plots, and saved as
#'   vector graphics with [save_life_cycle()] or [ggplot2::ggsave()].}
#'   \item{\code{"graphviz"}}{The original engine, which passes the matrix
#'   to \code{grViz} in \pkg{DiagrammeR} (see
#'   \url{https://rich-iannone.github.io/DiagrammeR/}) and returns an
#'   htmlwidget. Sub-matrices are summed into a single \code{matA}.}
#' }
#'
#' @param matA A matrix population model (i.e., a square projection matrix).
#'   Not needed if \code{matU} is supplied.
#' @param stages Optional vector of stage class labels. If missing, it first
#'   attempts to infer them from \code{dimnames(matA)} (or \code{matU}). If
#'   these are also \code{NULL}, then reverts to integers \code{1:ncol(A)}.
#' @param title Optional title for the plot. Defaults to \code{NULL}.
#' @param shape The shape to be used for the stages of the diagram. Any node
#'   shape accepted by \code{graphViz} is acceptable. Graphviz engine only.
#' @param fontsize Size of the font used in the diagram. Graphviz engine only.
#' @param nodefontsize Size of the font used in the node part of the diagram.
#'   Graphviz engine only.
#' @param edgecol Colour of the arrows in the diagram. Graphviz engine only.
#' @param node_order An optional numeric vector giving the position of each
#'   stage in the plot, so that stages are displayed in the order
#'   \code{stages[order(node_order)]}. Default is \code{NULL}, whereby the
#'   order is the same as \code{stages}, or row/column names, of the matrix.
#' @param matU,matF,matC Optional survival/growth (\code{matU}), sexual
#'   reproduction (\code{matF}) and clonal reproduction (\code{matC})
#'   matrices. Supply these instead of \code{matA} to colour transitions by
#'   type. \code{matF} and \code{matC} are optional.
#' @param engine Either \code{"ggplot"} (default) or \code{"graphviz"}.
#' @param ... Further options for the ggplot engine:
#'   \describe{
#'     \item{\code{digits}}{Decimal places for arrow labels (default 2).}
#'     \item{\code{show_values}}{Label arrows with their values (default
#'       \code{TRUE}).}
#'     \item{\code{scale_width}}{Scale arrow width by transition value,
#'       within each transition type (default \code{FALSE}).}
#'     \item{\code{colours}}{Named vector of colours for \code{"U"},
#'       \code{"F"}, \code{"C"}. Defaults are colour-blind safe.}
#'     \item{\code{type_labels}}{Named vector of legend labels for
#'       \code{"U"}, \code{"F"}, \code{"C"}.}
#'     \item{\code{node_fill}, \code{node_line}}{Fill and outline colour of
#'       stage nodes.}
#'     \item{\code{base_size}}{Base text size in pt (default 11). Nodes,
#'       arrows and the natural plot size scale with it.}
#'     \item{\code{base_family}}{Font family.}
#'   }
#'
#' @return With \code{engine = "ggplot"}, a \code{ggplot} object whose
#'   natural size (inches) is stored in \code{attr(, "lc_size")} and used by
#'   [save_life_cycle()]. With \code{engine = "graphviz"}, an object of class
#'   \code{grViz}.
#'
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#'
#' @family visualisation
#'
#' @examples
#' matA <- rbind(
#'   c(0.1, 0, 0, 0, 1.4),
#'   c(0.5, 0.2, 0, 0, 0),
#'   c(0, 0.3, 0.3, 0, 0),
#'   c(0, 0, 0.4, 0.4, 0.1),
#'   c(0, 0, 0, 0.1, 0.4)
#' )
#'
#' plot_life_cycle(matA)
#'
#' # Colour transitions by type by supplying the sub-matrices
#' data(mpm1)
#' plot_life_cycle(matU = mpm1$matU, matF = mpm1$matF)
#'
#' # Plot and label unknown transitions as NA
#' matA_na <- matA
#' matA_na[1, 5] <- NA
#' plot_life_cycle(matA_na)
#'
#' # Change the order of the nodes and give them names
#' plot_life_cycle(matA,
#'   stages = c("A", "B", "C", "D", "E"),
#'   node_order = 5:1
#' )
#'
#' # Larger text for slides, and arrow width scaled by value
#' plot_life_cycle(matU = mpm1$matU, matF = mpm1$matF,
#'                 base_size = 16, scale_width = TRUE)
#'
#' \dontrun{
#' # Save at natural size as vector graphics (PDF/SVG) or PNG
#' p <- plot_life_cycle(matU = mpm1$matU, matF = mpm1$matF)
#' save_life_cycle(p, "life_cycle.pdf")
#'
#' # The original graphviz engine
#' p1 <- plot_life_cycle(matA, engine = "graphviz")
#' p1 %>%
#'   DiagrammeRsvg::export_svg() %>%
#'   charToRaw() %>%
#'   rsvg::rsvg_png("my life cycle.png")
#' }
#'
#' @importFrom DiagrammeR grViz
#' @export plot_life_cycle
plot_life_cycle <- function(matA, stages, title = NULL, shape = "egg",
                            fontsize = 10, nodefontsize = 12,
                            edgecol = "grey", node_order = NULL,
                            matU = NULL, matF = NULL, matC = NULL,
                            engine = c("ggplot", "graphviz"), ...) {
  engine <- match.arg(engine)

  # Assemble matrices: either matA alone, or matU (+ matF, matC)
  if (!is.null(matU)) {
    if (!missing(matA)) {
      stop("Supply either 'matA' or 'matU' (with 'matF'/'matC'), not both",
           call. = FALSE)
    }
    mats <- list(U = matU, F = matF, C = matC)
  } else {
    if (missing(matA)) stop("Supply 'matA' or 'matU'", call. = FALSE)
    if (!is.null(matF) || !is.null(matC)) {
      stop("'matF' and 'matC' must be supplied with 'matU', not 'matA'",
           call. = FALSE)
    }
    mats <- list(U = matA)
  }
  mats <- mats[!vapply(mats, is.null, logical(1))]
  n <- nrow(mats[[1]])
  for (m in mats) {
    if (!is.matrix(m) || nrow(m) != ncol(m) || nrow(m) != n) {
      stop("All matrices must be square and of the same dimension",
           call. = FALSE)
    }
  }

  if (engine == "graphviz") {
    mat_sum <- Reduce(`+`, mats)
    if (missing(stages)) {
      return(.plot_life_cycle_graphviz(mat_sum, title = title, shape = shape,
                                       fontsize = fontsize,
                                       nodefontsize = nodefontsize,
                                       edgecol = edgecol,
                                       node_order = node_order))
    }
    return(.plot_life_cycle_graphviz(mat_sum, stages, title = title,
                                     shape = shape, fontsize = fontsize,
                                     nodefontsize = nodefontsize,
                                     edgecol = edgecol,
                                     node_order = node_order))
  }

  # Stage labels for the ggplot engine, mirroring the graphviz logic
  if (missing(stages)) {
    dn <- dimnames(mats[[1]])
    if (is.null(dn) || is.null(dn[[1]])) {
      stages <- seq_len(n)
    } else {
      stages <- dn[[1]]
      if (!identical(dn[[1]], dn[[2]])) {
        message(strwrap(
          prefix = " ", initial = "",
          "Dimension names of the matrix are not identical
      for rows and columns. Using row names."
        ))
      }
    }
  }
  if (length(stages) != n) {
    stop("The length of stages does not equal the dimension of the matrix",
         call. = FALSE)
  }
  if (!is.null(node_order) && length(node_order) != n) {
    stop("The length of node_order does not equal the dimension of the matrix",
         call. = FALSE)
  }

  .plot_life_cycle_gg(mats, stages = stages, node_order = node_order,
                      title = title, ...)
}


#' Life cycle diagram drawn with DiagrammeR/graphviz (internal engine)
#' @noRd
.plot_life_cycle_graphviz <- function(matA, stages, title = NULL,
                                      shape = "egg", fontsize = 10,
                                      nodefontsize = 12, edgecol = "grey",
                                      node_order = NULL) {
  # Identify stages
  if (missing(stages) && is.null(dimnames(matA))) {
    stages <- seq_len(ncol(matA))
  } else if (missing(stages) && !is.null(dimnames(matA))) {
    stages <- dimnames(matA)[[1]]

    if (!identical(
      dimnames(matA)[[1]],
      dimnames(matA)[[2]]
    )) {
      message(strwrap(
        prefix = " ", initial = "",
        "Dimension names of 'matA' are not identical
      for rows and columns. Using row names."
      ))
    }
  }

  if (length(stages) != nrow(matA)) {
    stop("The length of stages does not equal the dimension of matA")
  }

  if (!is.null(node_order)) {
    if (length(node_order) != nrow(matA)) {
      stop("The length of node_order does not equal the dimension of matA")
    }
  }
  # Construct a "from" -> "to" graph dataset (edges)
  graph <- expand.grid(to = stages, from = stages)
  graph$trans <- c(matA)
  graph$label <- ifelse(is.na(graph$trans), "NA", format(round(graph$trans, 3)))

  # Subset to include positive transitions and unknown transitions
  graph <- graph[is.na(graph$trans) | graph$trans > 0, ]

  # Create vector of node names (add semicolon for use by graphViz)
  if (is.null(node_order)) {
    nodes <- paste(paste0("'", stages, "'"), collapse = " -> ")
  } else {
    stages <- stages[order(node_order)]
    nodes <- paste(paste0("'", stages, "'"), collapse = " -> ")
  }
  nodes <- paste0(nodes, " [style=invis]")

  # Manipulate minimim length of edge to make the plot pretty (experimental!)
  graph$min_len <- (as.numeric(graph$to) - as.numeric(graph$from)) * 3
  # Create the edges argument for graphviz by pasting commands together
  edges <- paste0(
    "'", graph$from, "'", " -> ", "'", graph$to, "'",
    "[minlen=", graph$min_len,
    ",fontsize=", fontsize,
    ",color=", edgecol,
    ",xlabel=", paste("\"", graph$label),
    "\"]\n",
    collapse = ""
  )

  graph_label <- if (is.null(title)) {
    ""
  } else {
    paste0("  labelloc=\"t\";\n  label=\"", title, "\"\n")
  }

  # The graphviz argument, pasted together
  grViz(
    paste(
      "
digraph {
  {
    graph[overlap=false, margin=0, pad=0];
    rank=same;
    node [shape=", shape, ", fontsize=", nodefontsize, "];",
      nodes,
      "
  }",
      "
",
      edges,
      graph_label,
      "}"
    )
  )
}
