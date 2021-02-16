#' Calculate shape of survival over age
#'
#' Calculates a 'shape' value of survival lifespan inequality by comparing the 
#' area under a survival curve (over age) with the area under a constant 
#' survival function.
#'
#' @param surv Either 1) a numeric vector describing a survival curve (lx), or
#'   2) a \code{data.frame} / \code{list} with one column / element titled 'lx'
#'   describing a survival curve, optionally a column / element 'x' containing
#'   age classes (each element a number representing the age at the start of the
#'   class).
#'   
#'   If x is not supplied, the function will assume age classes starting at 0
#'   with time steps of 1 unit. If x begins at 0 then \code{lx[1]} should equal
#'   1. If x ends at maximum longevity, then \code{lx[which.max(x)]} should
#'   equal 0; however it is possible to supply partial survivorship curves.
#' @param xmin,xmax The minimum and maximum age respectively over which to
#'   evaluate shape. If not given, these default to \code{min(x)} and
#'   \code{max(x)} respectively.
#' @param trunc logical determining whether to truncate life tables or not 
#'   when any lx == 0. Usually this is the case only for the final value of lx.
#'   As the function calculates log(lx), these value(s) cannot be handled. 
#'   trunc == TRUE strips out the zero value(s). An alternative to this is to 
#'   transform the zeroes to something approximating zero (e.g. 1e-7).
#'
#' @return a shape value describing lifespan inequality by comparing the area
#'   under a survival (lx) curve over age with the area under a constant (type
#'   2) survival function. May take any real value between -0.5 and +0.5. A
#'   value of 0 indicates negligible ageing (neither generally increasing nor
#'   generally decreasing survival with age); negative values indicate negative
#'   senescence (generally increasing survival with age); positive values
#'   indicate senescence (generally decreasing survival with age). A value of
#'   +0.5 indicates that all individuals die at age of maximum longevity; a
#'   value of -0.5 indicates that (hypothetically) all individuals die at birth.
#' 
#' @author Iain Stott <iainmstott@@gmail.com>
#' 
#' @examples
#' # exponential decline in lx yields shape = 0
#' lx <- 0.7^(0:20)
#' shape_surv(lx)
#' 
#' @export shape_surv
shape_surv <- function(surv, xmin = NULL, xmax = NULL, trunc = FALSE) {
  if(class(surv) %in% "numeric") {
    lx <- surv
    x <- seq_along(lx) - 1
    if(lx[1] != 1) {
      stop("if x isn't given, lx must start with 1 as x[1] is assumed to be 0")
    }
  }
  if(class(surv) %in% c("list", "data.frame")) {
    if(!all(c("x", "lx") %in% names(surv))) {
      stop("'surv' doesn't contain both x and lx")
    }
    x <- surv$x
    lx <- surv$lx
    if(length(x) != length(lx)) {
      stop("x and lx must be the same length")
    }
    if((x[1] %in% 0) & !(lx[1] %in% 1)){
      stop("lx must start with 1 where x[1] is 0")
    }
  }
  if(!trunc) {
    if(any(lx %in% 0)) {
      stop(
"lx cannot be zero (we calculate the log). Consider trunc = TRUE, \n
or transforming zero values. See ?shape_surv for more details.\n")
    }
  }
  if(trunc) {
    x <- x[lx > 0]
    lx <- lx[lx > 0]
  }
  if(is.null(xmin)) xmin <- min(x)
  if(is.null(xmax)) xmax <- max(x)
  if(any(diff(x) <= 0)) stop("much as we'd like to reverse ageing, x must all be ascending")
  if(any(diff(lx) > 1e-7)) stop("please don't bring people back from the dead (check lx)")
  x_sub <- x[x >= xmin & x <= xmax]
  if(length(x_sub) <= 2) {
    stop("must have > 2 nonzero values of lx to calculate shape")
  }
  lx_sub <- lx[x >= xmin & x <= xmax]
  lx_log <- log(lx_sub)
  xStd <- (x_sub - xmin) / (xmax - xmin)
  lxmin <- lx_log[which.min(xStd)]
  lxmax <- lx_log[which.max(xStd)]
  lxStd <- (lx_log - lxmin) / (lxmax - lxmin)
  aucStd <- area_under_curve(xStd, lxStd)
  aucFlat <- 0.5
  shape <- aucFlat - aucStd
  shape
}
