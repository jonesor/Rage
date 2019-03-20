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
#'   In all cases where x ends at maximum longevity, \code{lx[which.max(x)]} 
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
#'
#' @export shape_rep
shape_rep <- function(rep, xmin = NULL, xmax = NULL, ...) {
  if(class(rep) %in% "numeric") {
    mx <- rep
    ltdim <- length(mx)
    x <- 0:(ltdim - 1)
    lt <- data.frame(x, mx)
  }
  if(class(rep) %in% "data.frame") {
    if(!all(c("x", "mx") %in% names(rep))) {
      stop("'rep' doesn't contain both x and mx")
    }
    lt <- rep[, c("x", "mx")]
    ltdim <- dim(lt)[1]
  }
  listtype <- "none"
  if(class(rep) %in% "list") {
    if(!all(c("x", "mx") %in% names(rep)) & !all(c("matU", "matF") %in% names(rep))) {
      stop("Please pass EITHER 'x' AND 'mx' OR 'matU' AND 'matF' to 'rep'")
    }
    if(all(c("x", "mx") %in% names(rep))) {
      if(any(c("matU", "matF") %in% names(rep))){
        stop("Please pass EITHER 'x' AND 'mx' OR 'matU' AND 'matF' to 'rep'")
      }
      if(length(unique(lengths(rep[c("x", "mx")]))) != 1) {
        stop("x and mx must be the same length")
      }
      lt <- as.data.frame(rep[c("x", "mx")])
      ltdim <- dim(lt)[1]
      listtype <- "xmx"
    }
    if(all(c("matU", "matF") %in% names(rep))) {
      if(any(c("x", "mx") %in% names(rep))){
        stop("Please pass EITHER 'x' AND 'mx' OR 'matU' AND 'matF' to 'rep'")
      }
      matU <- rep$matU
      matF <- rep$matF
      listtype <- "UF"
    }
  }
  if(class(rep) %in% c("list", "CompadreMat") & any(listtype %in% c("none", "UF"))){
    if(class(rep) %in% "CompadreMat") {
      matU <- Rcompadre::matU(rep)
      matF <- Rcompadre::matF(rep)
    }
    dots <- list(...)
    mLTargs <- c(list(matU = matU, matF = matF), dots[!names(dots) %in% "conv"])
    lt0 <- do.call("makeLifeTable", mLTargs)
    qC <- qsdConverge(matU, ...)
    mx <- lt0$mx[1:qC]
    mx[qC] <- 0
    ltdim <- qC
    x <- 0:(ltdim - 1)
    lt <- data.frame(x, mx)
  }
  if(is.null(xmin)) xmin <- min(which(lt$mx > 0))
  if(is.null(xmax)) xmax <- max(lt$x)
  if((xmax - xmin) <= 1) stop("xmax - xmin must be larger than 1")
  if(any(duplicated(lt$x))) stop("all x must be unique values")
  if(any(diff(lt$x) <= 0)) stop("x must all be ascending")
  lt$Bx <- c(0, cumsum(lt$mx[1:(ltdim - 1)]))
  if(any(diff(lt$x) < 0)) stop("no negative fecundity, thank you very much (check mx)")
  xStd <- (lt$x - xmin) / (xmax - xmin)
  Bxmin <- lt$Bx[which(lt$x %in% xmin)]
  Bxmax <- lt$Bx[which(lt$x %in% xmax)]
  BxStd <- (lt$Bx - Bxmin) / (Bxmax - Bxmin)
  aucStd <- .RageAUC(xStd, BxStd)
  aucFlat <- 0.5
  shape <- aucFlat - aucStd
  shape
}
