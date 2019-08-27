#' Calculate longevity from a matrix population model
#'
#' Calculate longevity (the age a which survivorship falls to some critical
#' proportion) from a matrix population model
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param start The index of the first stage at which the author considers the
#'   beginning of life. Defaults to 1. Alternately, a numeric vector giving the
#'   starting population vector (in which case \code{length(start)} must match
#'   \code{ncol(matU))}. See section \emph{Starting from multiple stages}.
#' @param xmax The maximum age to which survivorship will be calculated.
#'   Defaults to 1000.
#' @param lx_crit The critical proportion to calculate longevity with respect to
#'   (a value between 0 and 1). Defaults to 0.01.
#' 
#' @return Returns longevity, the integer age at which expected survivorship
#'   falls below \code{lx_crit}. If survivorship doesn't reach \code{lx_crit} by
#'   \code{xmax}, returns \code{NA} and prints a warning message.
#' 
#' @section Starting from multiple stages:
#' Rather than specifying argument \code{start} as a single stage class from
#' which all individuals start life, it may sometimes be desirable to allow for
#' multiple starting stage classes. For example, if we want to start our
#' calculation of longevity from reproductive maturity (i.e. first
#' reproduction), we should account for the possibility that there may be
#' multiple stage classes in which an individual could first reproduce.
#' 
#' To specify multiple starting stage classes, specify argument \code{start} as
#' the desired starting population vector (\strong{n1}), giving the proportion
#' of individuals starting in each stage class (the length of \code{start}
#' should match the number of columns in the relevant MPM).
#' 
#' See function \code{\link{mature_distrib}} for calculating the proportion of
#' individuals achieving reproductive maturity in each stage class.
#'  
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <hcaswell@@whoi.edu>
#' 
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'
#'   Morris, W. F., and D. F. Doak. (2003) Quantitative Conservation Biology:
#'   Theory and Practice of Population Viability Analysis. Sinauer Associates,
#'   Sunderland, Massachusetts, USA
#' 
#' @examples
#' data(mpm1)
#' 
#' longevity(mpm1$matU, start = 2)
#' longevity(mpm1$matU, start = 2, lx_crit = 0.05)
#' 
#' ### starting from first reproduction
#' repstages <- repro_stages(mpm1$matF)
#' n1 <- mature_distrib(mpm1$matU, start = 2, repro_stages = repstages)
#' 
#' longevity(mpm1$matU, start = n1)
#' 
#' @export longevity
longevity <- function(matU, start = 1L, xmax = 1000, lx_crit = 0.01) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU, start_vec = TRUE)
  if (lx_crit < 0 | lx_crit > 1) {
    stop("lx_crit must be a proportion between 0 and 1", call. = FALSE)
  }
  
  if (length(start) > 1) {
    n <- start
  } else {
    n <- rep(0.0, nrow(matU))
    n[start] <- 1.0
  }
  
  lx <- sum(n)
  t <- 0L
  
  while (lx > lx_crit & t < xmax) {
    n <- matU %*% n
    lx <- sum(n)
    t <- t + 1L
  }
  
  if (lx <= lx_crit) {
    longevity <- t
  } else {
    longevity <- NA_integer_
    warning("survivorship did not reach lx_crit by xmax: returning NA")
  }
  
  return(longevity)
}
