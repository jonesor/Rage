#' Calculate Keyfitz's entropy from a trajectory of age-specific survivorship
#' 
#' Calculate Keyfitz's entropy from a vector of age-specific survivorship (lx).
#' 
#' @section Warning:
#' Note that this function may produce unexpected results if used on partial 
#' survivorship trajectories. In addition, it is sensitive to the length of the 
#' survivorship vector. We direct users to the function `\code{\link{shape_surv}}` 
#' which is relatively robust to these issues.
#' 
#' @param lx Survivorship trajectory (a vector of monotonically-declining values
#'   in the interval [0,1]).
#' @param trapeze A logical argument indicating whether the composite trapezoid
#'   approximation should be used for approximating the definite integral.
#'   
#' @return Keyfitz's life table entropy.
#' 
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' 
#' @references Keyfitz, N. 1977. Applied Mathematical Demography. New York:
#'   Wiley.
#'   
#'   Demetrius, L., & Gundlach, V. M. 2014. Directionality theory and
#'   the entropic principle of natural selection. Entropy 16: 5428-5522.
#' 
#' @family {life history traits}
#' 
#' @examples
#' data(mpm1)
#' 
#' # derive lx trajectory, starting from stage 2
#' lx <- mpm_to_lx(mpm1$matU, start = 2)
#' 
#' # calculate Keyfitz' entropy
#' entropy_k(lx)
#' 
#' # use trapezoid approximation for definite integral
#' entropy_k(lx, trapeze = TRUE)
#' 
#' @export entropy_k
entropy_k <- function(lx, trapeze = FALSE) {

  # validate arguments
  if (any(lx < 0 | lx > 1)) {
    stop("All values of lx must be within the interval [0, 1]")
  }
  if (any(diff(lx) > 1e-7)) {
    stop("Values of lx must be monotonically declining")
  }
  
  # remove ages at/beyond which lx is 0 or NA
  lx <- lx[1:max(which(lx > 0))]
  
  # Calculate Keyfitz's entropy
  if (trapeze == TRUE) {
    x <- seq_along(lx)
    H <- -area_under_curve(x, lx * log(lx)) / area_under_curve(x, lx)
  } else {
    H <- -sum(lx * log(lx)) / sum(lx)
  }
  
  return(H) 
}
