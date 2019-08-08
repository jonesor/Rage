#' Calculate Keyfitz's entropy from a trajectory of age-specific survivorship
#' 
#' This function calculates Keyfitz's entropy from a vector of age-specific
#' survivorship (lx).
#' 
#' @param lx Survivorship trajectory (a vector of monotonically-declining values
#'   in the interval [0,1])
#' @param trapeze A logical argument indicating whether the composite trapezoid
#'   approximation should be used for approximating the definite integral.
#'   
#' @return Keyfitz's life table entropy.
#' 
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' 
#' @references Keyfitz, N. (1977) Applied Mathematical Demography. New York:
#'   Wiley.
#' 
#' @examples
#' lx <- 0.8^(0:20)
#' 
#' entropy_k(lx)
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
