#' Calculate Demetrius' entropy from a matrix population model
#'
#' Calculates Demetrius' entropy from a matrix population model, by first using
#' age-from-stage decomposition methods to estimate age-specific survivorship
#' (lx) and fecundity (mx).
#' 
#' @param matU A square matrix containing only survival-related transitions
#'   (i.e. progression, stasis, retrogression).
#' @param matR A square matrix containing only reproduction-related transitions
#'   (either sexual, clonal, or both; i.e. \code{matF}, \code{matC}, or
#'   \code{matF + matC}).
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param nSteps The maximum age to which age-specific survival (lx) and
#'   reproduction (mx) will be calculated. This allows the user to exclude ages
#'   after which mortality or fertility has plateaued (see function
#'   \code{\link{qsdConverge}} for more information). Defaults to 100.
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
#' dEntropy(matU, matF, nSteps = 10)
#' dEntropy(matU, matF, nSteps = 100)
#' 
#' @export dEntropy
dEntropy <- function(matU, matR, startLife = 1, nSteps = 100) {
  
  # validate arguments
  if (any(is.na(matU))) {
    stop("matU contains NAs")
  } 
  if (any(is.na(matR))) {
    stop("matR contains NAs")
  } 
  if (any(colSums(matU) > 1)) {
    warning("matU has at least one stage-specific survival value > 1")
  }

  # age-specific survivorship (lx) and fecundity (mx)
  lx <- ageSpecificSurv(matU, startLife, nSteps)
  mx <- ageSpecificRepro(matU, matR, startLife, nSteps)
  
  # calculate Demetrius' entropy
  H <- dEntropyCalc(lx, mx)
  
  return(H)
}


# Calculate Demetrius' entropy given lx and mx
dEntropyCalc <- function(lx, mx) {
  lxmx <- lx * mx

  # if lxmx == 0, log(lxmx) == -Inf; for entropy calc below, these -Inf can be
  #  converted to 0, because lim(x * log(x)) as x->0 is 0
  log_lxmx <- log(lxmx)
  log_lxmx[lxmx == 0] <- 0

  H <- -sum(lxmx * log_lxmx) / sum(lxmx)
  return(H)
}
