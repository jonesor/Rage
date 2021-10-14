#' Add stage names to matrices
#' 
#' Adds user-supplied or automatically-generated stage names to a matrix
#' population model (MPM).
#'
#' @param mat An MPM, either as a single matrix or list of matrices.
#' @param names A character vector specifying the name of each life stage, in
#'   order. If provided, \code{prefix} and \code{left_pad} arguments are ignored.
#' @param prefix A string to be pre-pended to the stage number when automatically
#'   naming stages. Defaults to \code{stage_}.
#' @param left_pad Logical, whether to pre-pend \code{0} to stage names such
#'   that all stage numbers have equal length, enabling lexicographic sorting.
#'   For example, stage \code{1} becomes \code{01} for matrices with 10-99
#'   stages, \code{001} for matrices with 100-999 stages, and so on. Defaults to
#'   \code{TRUE}.
#'   
#' @author William K. Petry <wpetry@@ncsu.edu>
#' 
#' @family transformation
#'
#' @return The input matrix or matrices with named rows and columns.
#' @export name_stages
#'
#' @examples
#' matU <- rbind(c(0.0, 0.0, 0.0),
#'               c(0.3, 0.1, 0.0),
#'               c(0.0, 0.5, 0.8))
#' # (semi)automated naming
#' name_stages(matU)
#' name_stages(matU, prefix = "s")
#' # custom stage names
#' name_stages(matU, names = c("small", "medium", "large"))
#' # overwrite existing stage names
#' data(mpm1)
#' name_stages(mpm1)

name_stages <- function(mat, names = NULL, prefix = "stage_", left_pad = TRUE) {
  # check argument inputs
  if (!is.matrix(mat) && !is.list(mat)) {
    stop("Argument `mat` must be either a matrix or list of matrices.\n")
  }
  if (is.matrix(mat)) {
    if (length(unique(dim(mat))) != 1L) {
      stop("When `mat` is supplied as a matrix, it must be square.\n")
    } else {
      mdim <- nrow(mat)
    }
  } else {
    if (length(unique(unlist(lapply(mat, dim), use.names = FALSE))) != 1L) {
      stop("Each matrix in `mat` must be square with the same dimensions.\n")
    } else {
      mdim <- nrow(mat[[1]])
    }
  }
  if (is.null(names) && !is.character(prefix)) {
    stop("Either stage `names` or a naming `prefix` must be supplied.\n")
  }
  if (!is.null(names) && is.character(prefix)){
    warning("Naming `prefix` ignored, using stage `names` instead.\n")
  }
  # construct stage names
  if (is.null(names)) {
    if (isTRUE(left_pad)) {
      fmt <- paste0("%0", nchar(mdim), "d")
      names <- paste0(prefix, sprintf(fmt, 1:mdim))
    } else if (isFALSE(left_pad)) {
      names <- 1:mdim
    } else {
      stop("Argument `left_pad` must be of type logical (TRUE/FALSE).\n")
    }
  }
  if (mdim != length(names)) {
    stop("Incorrect number of stage names supplied: ", length(names),
         " names supplied for ", mdim, " life stages.\n")
  }
  # warn if overwriting existing stage names
  if (is.list(mat) && !is.null(unlist(lapply(mat, dimnames), use.names = FALSE)) |
      is.matrix(mat) && !is.null(unlist(dimnames(mat), use.names = FALSE))) {
    warning("Existing stage names have been overwritten!\n")
  }
  # add stage names to matrix/matrices
  if (is.matrix(mat)) {
    dimnames(mat) <- list(names, names)
  } else {
    mat <- lapply(mat, `dimnames<-`, list(names, names))
  }
  return(mat)
}
