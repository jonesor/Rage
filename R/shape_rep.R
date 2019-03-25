#' Calculate shape of reproduction over age
#'
#' Calculates a 'shape' value of distribution of reproduction over age by comparing
#' the area under a cumulative reproduction curve (over age) with the area 
#' under a cumulative function describing constant reproduction.
#'
#' @param rep Either 1) a numeric vector describing reproduction over age (mx); 2) a 
#'   \code{data.frame} / \code{list} with one column / element titled 'mx' 
#'   describing a reproduction over age, optionally a column / element 'x' containing 
#'   age classes (each element a number representing the age at the start of the
#'   class); 3) a list containing two elements: 'matU', a U matrix (the survival 
#'   component of a matrix population model, i.e. a square projection matrix reflecting 
#'   survival-related transitions, e.g. progression, stasis, and retrogression) and
#'   'matF', and F matrix(the reproduction component of a matrix projection model, 
#'   i.e. a square projection matrix of the same dimension as matU reflecting 
#'   those transitions describing sexual reproduction); 4) a \code{CompadreMat} 
#'   object (\code{RCompadre-package})containing a matrix population model in 
#'   the format described in the \code{CompadreMat} class.
#'   In the case of 1 and 2 where x is not supplied, the function will assume
#'   age classes starting at 0 with steps of 1 year.
#'   In the case of 3 and 4, an age-based reproduction schedule will be generated from 
#'   a stage-based matrix using the \code{makeLifeTable} function of 
#'   \code{RCompadre-package}. 
#'   In all cases where x ends at maximum longevity, \code{mx[which.max(x)]} 
#'   should equal 0, however it is possible to supply partial reproduction schedules.
#' @param xmin,xmax The minimum and maximum age repectively over which to evaluate
#'   shape. If not given, these default to \code{min(x)} and \code{max(x)} 
#'   respectively.
#' @param ... when \code{rep} is either U and F matrices or \code{CompadreMat} object,
#'   \code{...} specifies further paramters to pass to 
#'   \code{\link{makeLifeTable}} and to pass to \code{link{qsdConverge}}. 
#'   Can take \code{nsteps}, \code{startLife} and \code{conv}; see 
#'   \code{\link{makeLifeTable}} and \code{link{qsdConverge}}, as it may be
#'   important to adjust these for your model in order to generate a 
#'   meaningful life table.
#'
#' @return a shape value describing symmetry of reproduction over age by comparing 
#'   the area under a cumulative reproduction curve over age with the area under 
#'   constant reproduction. May take any real value between -0.5 and +0.5. A value
#'   of 0 indicates negligible aging (neither generally increasing nor generally
#'   decreasing reproduction with age); negative values indicate
#'   senescence (generally decreasing reproduction with age); positive values indicate 
#'   negative senescence (generally increasing reproduction with age). A value 
#'   of -0.5 indicates that (hypothetically) all individuals are born to 
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
shape_rep <- function(rep, xmin = NULL, xmax = NULL, ...) {
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
  if(is.null(xmin)) xmin <- min(which(mx > 0))
  if(is.null(xmax)) xmax <- max(x)
  if(any(duplicated(x))) stop("all x must be unique values")
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
  Bxmin <- Bx[which.min(xStd)]
  Bxmax <- Bx[which.max(xStd)]
  BxStd <- (Bx_sub - Bxmin) / (Bxmax - Bxmin)
  aucStd <- area_under_curve(xStd, BxStd)
  aucFlat <- 0.5
  shape <- aucStd - aucFlat
  shape
}
