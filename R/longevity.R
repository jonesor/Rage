#' Calculate longevity from a matrix population model
#'
#' Calculate longevity (the age a which survivorship falls to some critical
#' proportion) from a matrix population model
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param startLife Index of the first stage at which the author considers the
#'   beginning of life. Defaults to 1.
#' @param lxCrit The critical proportion to calculate longevity with respect to
#'   (a value between 0 and 1). Defaults to 0.01.
#' @param maxAge The maximum age to which survivorship will be calculated.
#'   Defaults to 1000.
#' @return Returns longevity, the integer age at which expected survivorship
#'   falls below \code{lxCrit}. If survivorship doesn't reach \code{lxCrit} by
#'   \code{maxAge}, returns \code{NA} and prints a warning message.
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <hcaswell@@whoi.edu>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'
#'   Morris, W. F., and D. F. Doak. (2003) Quantitative Conservation Biology:
#'   Theory and Practice of Population Viability Analysis. Sinauer Associates,
#'   Sunderland, Massachusetts, USA
#' @examples
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#'
#' longevity(matU)
#' longevity(matU, lxCrit = 0.05)
#' longevity(matU, startLife = 3, lxCrit = 0.05, maxAge = 50)
#' 
#' @export longevity
longevity <- function(matU, startLife = 1, lxCrit = 0.01, maxAge = 1000) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(startLife, matU)
  if (lxCrit < 0 | lxCrit > 1) {
    stop("lxCrit must be a proportion between 0 and 1", call. = FALSE)
  }
  
  matUtemp <- matU
  lx <- 1.0
  t <- 0L
  
  while (lx > lxCrit & t < maxAge) {
    lx <- sum(matUtemp[,startLife])
    matUtemp <- matUtemp %*% matU
    t <- t + 1L
  }
  
  if (lx <= lxCrit) {
    longevity <- t
  } else {
    longevity <- NA_integer_
    warning("survivorship did not reach lxCrit by maxAge: returning NA")
  }
  
  return(longevity)
}
