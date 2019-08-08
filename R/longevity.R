#' Calculate longevity from a matrix population model
#'
#' Calculate longevity (the age a which survivorship falls to some critical
#' proportion) from a matrix population model
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param start Index of the first stage at which the author considers the
#'   beginning of life. Defaults to 1.
#' @param xmax The maximum age to which survivorship will be calculated.
#'   Defaults to 1000.
#' @param lx_crit The critical proportion to calculate longevity with respect to
#'   (a value between 0 and 1). Defaults to 0.01.
#' 
#' @return Returns longevity, the integer age at which expected survivorship
#'   falls below \code{lx_crit}. If survivorship doesn't reach \code{lx_crit} by
#'   \code{xmax}, returns \code{NA} and prints a warning message.
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
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#'
#' longevity(matU)
#' longevity(matU, lx_crit = 0.05)
#' longevity(matU, start = 3, xmax = 50, lx_crit = 0.05)
#' 
#' @export longevity
longevity <- function(matU, start = 1, xmax = 1000, lx_crit = 0.01) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU)
  if (lx_crit < 0 | lx_crit > 1) {
    stop("lx_crit must be a proportion between 0 and 1", call. = FALSE)
  }
  
  tempU <- matU
  lx <- 1.0
  t <- 0L
  
  while (lx > lx_crit & t < xmax) {
    lx <- sum(tempU[,start])
    tempU <- tempU %*% matU
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
