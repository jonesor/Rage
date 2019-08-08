#' Calculate Demetrius' entropy from trajectories of age-specific survivorship
#' and fecundity
#'
#' This function calculates Demetrius' entropy from vectors of age-specific
#' survivorship (lx) and fecundity (mx).
#' 
#' @param lx Age-specific survivorship trajectory (a vector of
#'   monotonically-declining values in the interval [0,1])
#' @param mx Age-specific fecundity trajectory (a vector of non-negative values)
#'   
#' @return Demetrius' entropy.
#' 
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @references Demetrius, L. (1978) Adaptive value, entropy and survivorship
#'   curves. Nature 275, 213-214. doi:10.1038/275213a0
#'   
#' @examples
#' lx <- c(1.00, 0.85, 0.70, 0.65, 0.55, 0.50, 0.45, 0.40, 0.35)
#' mx <- c(0.00, 0.00, 1.10, 1.50, 1.60, 1.70, 1.50, 1.20, 0.70)
#' 
#' entropy_d(lx, mx)
#' 
#' @export entropy_d
entropy_d <- function(lx, mx) {
  
  # validate arguments
  if (any(lx < 0 | lx > 1)) {
    stop("All values of lx must be within the interval [0, 1]")
  }
  if (any(diff(lx) > 1e-7)) {
    stop("Values of lx must be monotonically declining")
  }
  if (any(mx < 0)) {
    stop("All values of mx must be >= 0")
  }
  
  # calculate Demetrius' entropy
  lxmx <- lx * mx
  # if lxmx == 0, log(lxmx) == -Inf; for entropy calc below, these -Inf can be
  #  converted to 0, because lim(x * log(x)) as x->0 is 0
  log_lxmx <- log(lxmx)
  log_lxmx[lxmx == 0] <- 0
  
  H <- -sum(lxmx * log_lxmx) / sum(lxmx)
  
  return(H)
}
