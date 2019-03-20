#' Calculate shape of survival over age
#'
#' Calculates a 'shape' value of survival lifespan inequality by comparing the 
#' area under a survival curve (over age) with the area under a constant 
#' survival function.
#'
#' @param surv Either 1) a numeric vector describing a survival curve (lx); 2) a 
#'   \code{data.frame} / \code{list} with one column / element titled 'lx' 
#'   describing a survival curve, optionally a column / element 'x' containing 
#'   age classes (each element a number representing the age at the start of the
#'   class); 3) a 'U' matrix (the survival component of a matrix 
#'   population model, i.e. a square projection matrix reflecting 
#'   survival-related transitions, e.g. progression, stasis, and retrogression); 
#'   4) a \code{CompadreM} object (\code{RCompadre-package})containing a matrix 
#'   population model in the format described in the \code{CompadreM} class.
#'   In the case of 1 and 2 where x is not supplied, the function will assume
#'   age classes starting at 0 with steps of 1 year.
#'   In the case of 3 and 4, an age-based survival curve will be generated from 
#'   a stage-based matrix using \code{mpm_to_lx}. 
#'   In all cases where x begins at 0 and ends at maximum longevity,
#'   \code{lx[1]} should equal 1 and \code{lx[which.max(x)]} should equal 0,
#'   however it is possible to supply partial survival curves.
#' @param xmin,xmax The minimum and maximum age respectively over which to
#'   evaluate shape. If not given, these default to \code{min(x)} and
#'   \code{max(x)} respectively.
#' @param ... when \code{surv} is either a U matrix or \code{CompadreM} object,
#'   \code{...} specifies further parameters to pass to \code{\link{mpm_to_lx}}
#'   and \code{link{qsdConverge}}. Can take \code{nSteps}, \code{startLife} and
#'   \code{conv}; see \code{\link{mpm_to_lx}} and \code{link{qsdConverge}}, as
#'   it may be important to adjust these for your model in order to generate a
#'   meaningful life table.
#'
#' @return a shape value describing lifespan inequality by comparing the area
#'   under a survival (lx) curve over age with the area under a constant (type
#'   2) survival function. May take any real value between -0.5 and +0.5. A
#'   value of 0 indicates negligible aging (neither generally increasing nor
#'   generally decreasing survival with age); negative values indicate
#'   senescence (generally decreasing survival with age); positive values
#'   indicate negative senescence (generally increasing survival with age). A
#'   value of -0.5 indicates that all individuals die at age of maximum
#'   longevity; a value of +0.5 indicates that (hypothetically) all individuals
#'   die at birth.
#' 
#' @author Iain Stott <iainmstott@@gmail.com>
#' 
#' @examples
#' # exponential decline in lx yields shape = 0
#' lx <- 0.7^(0:20)
#' shape_surv(lx)
#' 
#' @importMethodsFrom Rcompadre matU
#' @export shape_surv
shape_surv <- function(surv, xmin = NULL, xmax = NULL, ...) {
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
  if(class(surv) %in% c("matrix", "CompadreMat")){
    if(class(surv) %in% "CompadreMat") {
      matU <- matU(surv)
    } else {
      matU <- surv
    }
    dots <- list(...)
    mLTargs <- c(list(matU = matU), dots[!names(dots) %in% "conv"])
    lx <- do.call("mpm_to_lx", mLTargs)
    qC <- qsdConverge(matU, ...)
    if (!(is.na(qC) | qC > length(lx))) lx <- lx[1:qC]
    x <- seq_along(lx) - 1
    if(lx[1] != 1) {
      stop("error in mpm_to_lx: lx[1] != 1")
    }
  }
  x <- x[lx > 0]
  lx <- lx[lx > 0]
  if(is.null(xmin)) xmin <- min(x)
  if(is.null(xmax)) xmax <- max(x)
  if(any(duplicated(x))) stop("all x must be unique values")
  if(any(diff(x) <= 0)) stop("x must all be ascending")
  if(any(diff(lx) > 1e-7)) {
    stop("please don't bring people back from the dead (check lx)")
  }
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
