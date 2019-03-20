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
#'   a stage-based matrix using the \code{makeLifeTable} function of 
#'   \code{RCompadre-package}. 
#'   In all cases where x begins at 0 and ends at maximum longevity, \code{lx[1]} 
#'   should equal 1 and \code{lx[which.max(x)]} should equal 0, however it is 
#'   possible to supply partial survival curves.
#' @param xmin,xmax The minimum and maximum age repectively over which to evaluate
#'   shape. If not given, these default to \code{min(x)} and \code{max(x)} 
#'   respectively.
#' @param ... when \code{surv} is either a U matrix or \code{CompadreM} object,
#'   \code{...} specifies further paramters to pass to 
#'   \code{\link{makeLifeTable}} and to pass to \code{link{qsdConverge}}. 
#'   Can take \code{nsteps}, \code{startLife} and \code{conv}; see 
#'   \code{\link{makeLifeTable}} and \code{link{qsdConverge}}, as it may be
#'   important to adjust these for your model in order to generate a 
#'   meaningful life table.
#'
#' @return a shape value describing lifespan inequality by comparing 
#'   the area under a survival (lx) curve over age with the area under a constant 
#'   (type 2) survival function. May take any real value between -0.5 and +0.5. A value
#'   of 0 indicates negligible aging (neither generally increasing nor generally
#'   decreasing survival with age); negative values indicate
#'   senescence (generally decreasing survival with age); positive values indicate negative senescence 
#'   (generally increasing survival with age). A value of -0.5 indicates that 
#'   all individuals die at age of maximum longevity; a value of +0.5 
#'   indicates that (hypothetically) all individuals die at birth.
#' 
#' @author Iain Stott <iainmstott@@gmail.com>
#' 
#' @examples
#'
#' @export shape_surv
shape_surv <- function(surv, xmin = NULL, xmax = NULL, ...) {
  if(class(surv) %in% "numeric") {
    lx <- surv
    ltdim <- length(lx)
    x <- 0:(ltdim - 1)
    lt <- data.frame(x, lx)
    if(!(lt$lx[1] %in% 1)){
      stop("if x isn't given, lx must start with 1 as x[1] is assumed to be 0")
    }
  }
  if(class(surv) %in% "data.frame") {
    if(!all(c("x", "lx") %in% names(surv))) {
      stop("'surv' doesn't contain both x and lx")
    }
    lt <- surv[, c("x", "lx")]
    ltdim <- dim(lt)[1]
    if((lt$x[1] %in% 0) & !(lt$lx[1] %in% 1)){
      stop("lx must start with 1 where x[1] is 0")
    }
  }
  if(class(surv) %in% "list") {
    if(!all(c("x", "lx") %in% names(surv))) {
      stop("'surv' doesn't contain both x and lx")
    }
    if(length(unique(lengths(surv[c("x", "lx")]))) != 1) {
      stop("x and lx must be the same length")
    }
    lt <- as.data.frame(surv[c("x", "lx")])
    ltdim <- dim(lt)[1]
    if((lt$x[1] %in% 0) & !(lt$lx[1] %in% 1)){
      stop("lx must start with 1 where x[1] is 0")
    }
  }
  if(class(surv) %in% c("matrix", "CompadreM")){
    if(class(surv) %in% "CompadreM") {
      matU <- RCompadre::matU(surv)
    } else {
      matU <- surv
    }
    dots <- list(...)
    mLTargs <- c(list(matU = matU), dots[!names(dots) %in% "conv"])
    lt0 <- do.call("makeLifeTable", mLTargs)
    qC <- qsdConverge(matU, ...)
    lx <- lt0$lx[1:qC]
    lx[qC] <- 0
    ltdim <- qC
    x <- 0:(ltdim - 1)
    
    lt <- data.frame(x, lx)
    if(!all(lt$lx[1] %in% 1)) {
      stop("error in makeLifeTable: lx[1] != 1")
    }
  }
  if(is.null(xmin)) xmin <- min(x)
  if(is.null(xmax)) xmax <- max(x)
  if((xmax - xmin) <= 1) stop("xmax - xmin must be larger than 1")
  if(any(duplicated(x))) stop("all x must be unique values")
  if(any(diff(x) <= 0)) stop("x must all be ascending")
  if(any(diff(lx) > 0)) stop("please don't bring people back from the dead (check lx)")
  xStd <- (x - xmin) / (xmax - xmin)
  lxmin <- lx[which(x %in% xmin)]
  lxmax <- lx[which(x %in% xmax)]
  lxStd <- (lx - lxmin) / (lxmax - lxmin)
  aucStd <- .RageAUC(xStd, lxStd)
  aucFlat <- 0.5
  shape <- aucFlat - aucStd
  shape
}

.RageAUC <- function(x, y, a = NULL, b = NULL) {
  if(is.null(a)){
    a <- which(x %in% min(x))
  } else {
    a <- which(x %in% a)
  }
  if(is.null(b)){
    b <- which(x %in% max(x))
  } else {
    b <- which(x %in% b)
  }
  if(any(diff(x) <=0)) stop("AUC: x should be ascending")
  if(a > b) stop("AUC: b should be greater than a")
  polys <- numeric()
  for(age in a:(b - 1)){
    polys[age] <- (x[age + 1] - x[age]) * min(c(y[age], y[age + 1])) + 
      0.5 * (x[age + 1] - x[age]) * (max(c(y[age], y[age + 1])) - min(c(y[age], y[age + 1])) )
  } 
  AUCab <- sum(polys[a:(b - 1)])
  AUCab
}
