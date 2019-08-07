#' Calculate Demetrius' entropy from a matrix population model
#'
#' Calculates Demetrius' entropy from a matrix population model, by first using
#' age-from-stage decomposition methods to estimate age-specific survivorship
#' (lx) and fecundity (mx).
#' 
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param nSteps The maximum age to which age-specific survival (lx) and
#'   reproduction (mx) will be calculated. This allows the user to exclude ages
#'   after which mortality or fertility has plateaued (see function
#'   \code{\link{qsd_converge}} for more information). Defaults to 100.
#' @return Demetrius' entropy.
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' @references Demetrius, L. (1978) Adaptive value, entropy and survivorship
#'   curves. Nature 275, 213-214. doi:10.1038/275213a0
#' @examples
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#' 
#' matF <- rbind(c(  0,   0, 1.1, 1.6),
#'               c(  0,   0, 0.8, 0.4),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#' 
#' entropy_d(matU, matF, nSteps = 10)
#' entropy_d(matU, matF, nSteps = 100)
#' 
#' @export entropy_d
entropy_d <- function(matU, matR, startLife = 1, nSteps = 100) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidMat(matR)
  checkValidStartLife(startLife, matU)

  # age-specific survivorship (lx) and fecundity (mx)
  lx <- mpm_to_lx(matU, startLife, nSteps)
  mx <- mpm_to_mx(matU, matR, startLife, nSteps)
  
  # calculate Demetrius' entropy
  H <- entropy_calc(lx, mx)
  
  return(H)
}


# Calculate Demetrius' entropy given lx and mx
entropy_calc <- function(lx, mx) {
  lxmx <- lx * mx

  # if lxmx == 0, log(lxmx) == -Inf; for entropy calc below, these -Inf can be
  #  converted to 0, because lim(x * log(x)) as x->0 is 0
  log_lxmx <- log(lxmx)
  log_lxmx[lxmx == 0] <- 0

  H <- -sum(lxmx * log_lxmx) / sum(lxmx)
  return(H)
}
