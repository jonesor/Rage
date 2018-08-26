#' Calculate longevity from a matrix population model
#'
#' Calculate longevity (the age a which survivorship falls to some critical
#' proportion) from a matrix population model
#'
#' @param matU A matrix containing only survival-dependent processes (growth,
#'   stasis, shrinkage).
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
#'   Morris, W. F., and D. F. Doak. 2003. Quantitative Conservation Biology:
#'   Theory and Practice of Population Viability Analysis. Sinauer Associates,
#'   Sunderland, Massachusetts, USA
#' @examples
#' matU <- rbind(c(0.0, 0.0, 0.0, 0.0),
#'               c(0.5, 0.0, 0.0, 0.0),
#'               c(0.0, 0.3, 0.0, 0.0),
#'               c(0.0, 0.0, 0.1, 0.1))
#'
#' longevity(matU)
#' longevity(matU, lxCrit = 0.05)
#' longevity(matU, startLife = 3, lxCrit = 0.05, maxAge = 50)
#' @export longevity
longevity <- function(matU, startLife = 1, lxCrit = 0.01, maxAge = 1000) {

  # validate inputs
  if (lxCrit < 0 | lxCrit > 1) {
    stop('lxCrit must be a proportion between 0 and 1')
  }
  if (any(is.na(matU))) {
    stop('matU contains missing values')
  }
  if (all(matU == 0)) {
    warning('all elements of matU are zero')
    return(1)
  }

  # calculate survivorship to maxAge
  lx <- ageSpecificSurv(matU, startLife, N = maxAge)
  
  if (any(lx <= lxCrit)) {
    longevity <- min(which(lx <= lxCrit)) - 1  # -1 b/c lx starts at age 0
  } else {
    longevity <- NA_real_
    warning('survivorship did not fall to lxCrit by maxAge: returning NA')
  }
  
	return(longevity)
}
