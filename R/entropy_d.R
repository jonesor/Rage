#' Calculate Demetrius' entropy from trajectories of age-specific survivorship
#' and fecundity
#'
#' This function calculates Demetrius' entropy from vectors of age-specific
#' survivorship (lx) and fecundity (mx).
#' 
#' #' @section Warning:
#' Note that this function may produce unexpected results if used on partial
#' survivorship and fecundity trajectories. In addition, it is sensitive to the
#' length of the these vectors. We direct users to the functions
#' `\code{\link{shape_surv}}` and `\code{\link{shape_rep}}` which are relatively
#' robust to these issues.
#' 
#' @param lx Age-specific survivorship trajectory (a vector of
#'   monotonically-declining values in the interval [0,1]).
#' @param mx Age-specific fecundity trajectory (a vector of non-negative values).
#'   
#' @return Demetrius' entropy.
#' 
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @family {life history traits}
#' 
#' @references Demetrius, L. 1978 Adaptive value, entropy and survivorship
#'   curves. Nature, 275, 213-214. <doi:10.1038/275213a0>
#'   
#' @examples
#' data(mpm1)
#' 
#' # derive trajectories of lx and mx, starting from stage 2
#' lx <- mpm_to_lx(mpm1$matU, start = 2)
#' mx <- mpm_to_mx(mpm1$matU, mpm1$matF, start = 2)
#' 
#' # calculate Demetrius' entropy
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
