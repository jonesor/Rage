#' Calculate shape of reproduction over age
#'
#' Calculates a 'shape' value of distribution of reproduction over age by
#' comparing the area under a cumulative reproduction curve (over age) with the
#' area under a cumulative function describing constant reproduction.
#'
#' @param rep Either 1) a numeric vector describing reproduction over age (mx),
#'   or 2) a \code{data.frame} / \code{list} with one column / element titled
#'   'mx' describing a reproduction over age, optionally a column / element 'x'
#'   containing age classes (each element a number representing the age at the
#'   start of the class).
#'   
#'   If x is not supplied, the function will assume age classes starting at 0
#'   with time steps of unit. If x ends at maximum longevity,
#'   \code{mx[which.max(x)]} should equal 0; however it is possible to supply
#'   partial reproduction schedules.
#' @param xmin,xmax The minimum and maximum age repectively over which to
#'   evaluate shape. If not given, these default to \code{min(x)} and
#'   \code{max(x)} respectively.
#'
#' @return a shape value describing symmetry of reproduction over age by
#'   comparing the area under a cumulative reproduction curve over age with the
#'   area under constant reproduction. May take any real value between -0.5 and
#'   +0.5. A value of 0 indicates negligible aging (neither generally increasing
#'   nor generally decreasing reproduction with age); negative values indicate
#'   senescence (generally decreasing reproduction with age); positive values
#'   indicate negative senescence (generally increasing reproduction with age).
#'   A value of -0.5 indicates that (hypothetically) all individuals are born to
#'   individuals of age 0; a value of +0.5 indicates that all individuals are
#'   born at the age of maximum longevity.
#' 
#' @author Iain Stott <iainmstott@@gmail.com>
#' 
#' @examples
#' mx <- c(0, 0, 0.3, 0.4, 0.5, 0.6)
#' shape_rep(mx)
#'
#' @export shape_rep
shape_rep <- function(rep, xmin = NULL, xmax = NULL) {
  if(class(rep) %in% "numeric") {
    mx <- rep
    x <- seq_along(mx) - 1
  }
  if(class(rep) %in% c("list", "data.frame")) {
    if(!all(c("x", "mx") %in% names(rep))) {
      stop("'rep' doesn't contain both x and mx")
    }
    x <- rep$x
    mx <- rep$mx
    if(length(x) != length(mx)) {
      stop("x and mx must be the same length")
    }
  }
  if(is.null(xmin)) xmin <- x[min(which(mx > 0))]
  if(is.null(xmax)) xmax <- max(x)
  if(any(diff(x) <= 0)) stop("much as we'd like to reverse aging, x must all be ascending")
  if(any(mx < 0)) stop("You appear to have minus-babies (check mx)")
  x_sub <- x[x >= xmin & x <= xmax]
  if(length(x_sub) <= 2) {
    stop("must have > 2 nonzero values of mx to calculate shape")
  }
  ltdim <- length(x)
  Bx <- c(0, cumsum(mx[1:(ltdim - 1)]))
  Bx_sub <- Bx[x >= xmin & x <= xmax]
  xStd <- (x_sub - xmin) / (xmax - xmin)
  Bxmin <- Bx_sub[which.min(xStd)]
  Bxmax <- Bx_sub[which.max(xStd)]
  BxStd <- (Bx_sub - Bxmin) / (Bxmax - Bxmin)
  aucStd <- area_under_curve(xStd, BxStd)
  aucFlat <- 0.5
  shape <- aucStd - aucFlat
  shape
}
