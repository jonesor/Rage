#' Plot a life cycle diagram from a matrix population model
#'
#' Plots the life cycle diagram illustrated by a matrix population model. This
#' function processes the matrix model and passes the information to the
#' graphViz function in DiagrammeR. See
#' \url{https://rich-iannone.github.io/DiagrammeR/}.
#'
#' @param matA A matrix population model (i.e., a square projection matrix)
#' @param stages Optional vector of stage class labels. If missing, it first
#'   attempts to infer them from \code{dimnames(matA)}. If these are also
#'   \code{NULL}, then reverts to integers \code{1:ncol(A)}.
#' @param title Optional title for the plot. Defaults to \code{NULL}.
#' @param shape The shape to be used for the stages of the diagram. Any node
#'   shape accepted by \code{graphViz} is acceptable.
#' @param fontsize Size of the font used in the diagram.
#' @param nodefontsize Size of the font used in the node part of the diagram.
#' @param edgecol Colour of the arrows in the diagram.
#' @param node_order An optional numeric vector giving the order that the nodes
#'   should be presented in the plot. Default is `NULL` whereby the order is the
#'   same as `stages`, or row/column names, of the matrix.
#'
#' @return An object of class \code{grViz} representing the life cycle diagram
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
#' # One could save the diagram as a PNG file using a combination of `export_svg`
#' # (from the `DiagrammeRsvg` package) and `rsvg_png` (from the `rsvg` package)
#' # like this:
#' \dontrun{
#' p1 <- plot_life_cycle(matA)
#' p1 %>%
#'   DiagrammeRsvg::export_svg %>%
#'   charToRaw() %>%
#'   rsvg::rsvg_png("my life cycle.png")
#' }
#'
#' # Change the order of the nodes and give them names
#' plot_life_cycle(matA,
#'   stages = c("A", "B", "C", "D", "E"),
#'   node_order = 5:1
#' )
#'
#' @importFrom DiagrammeR grViz
#' @export plot_life_cycle
plot_life_cycle <- function(matA, stages, title = NULL, shape = "egg",
                            fontsize = 10, nodefontsize = 12,
                            edgecol = "grey", node_order = NULL) {
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
  graph$trans <- round(c(matA), 3)

  # Subset to only include those where the trans > 0
  graph <- graph[graph$trans > 0, ]

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
  edges <- paste0("'", graph$from, "'", " -> ", "'", graph$to, "'",
    "[minlen=", graph$min_len,
    ",fontsize=", fontsize,
    ",color=", edgecol,
    ",xlabel=", paste("\"", graph$trans),
    "\"]\n",
    collapse = ""
  )

  # The graphviz argument, pasted together
  grViz(
    paste(
      "
digraph {
  {
    graph[overlap=false];
    rank=same;
    node [shape=", shape, ", fontsize=", nodefontsize, "];",
      nodes, "
  }",
      "ordering=out
  x [style=invis]
  x -> {", nodes, "} [style=invis]", edges,
      "labelloc=\"t\";
  label=\"", title, "\"
}"
    )
  )
}
