#' Calculate Demetrius' entropy
#' 
#' This function calculates Demetrius' entropy from a matrix population model,
#' by first using age-from-stage decomposition methods to estimate age-specific
#' survivorship (lx), fecundity (mx), and/or clonality (cx).
#' 
#' @param matU A square matrix containing only survival-related transitions
#'   (i.e. progression, stasis, retrogression).
#' @param matF A square matrix containing only sexual reproduction-related
#'   transitions.
#' @param matC A square matrix containing only clonal reproduction-related
#'   transitions.
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param nSteps The age-cutoff for the decomposition of age-specific survival
#'   (lx), sexual reproduction (mx), and/or clonal reproduction (cx).
#'   This allows the user to exclude ages after which mortality or fertility has
#'   plateaued (see function \code{qsdConverge} for more information). Defaults
#'   to 100.
#' @return A named list containing estimates of Demetrius' entropy with respect
#'   to sexual reproduction ('Fec', returned if \code{matF} is provided), clonal
#'   reproduction ('Clo', returned if \code{matC} is provided), and the sum of
#'   sexual and clonal reproduction ('FecClo', returned if both \code{matF} and
#'   \code{matC} are provided).
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @references Demetrius, L. (1978) Adaptive value, entropy and survivorship
#'   curves. Nature 275, 213-214. doi:10.1038/275213a0
#' @examples
#' matU <- rbind(c(0.0, 0.0, 0.0, 0.0),
#'               c(0.5, 0.0, 0.0, 0.0),
#'               c(0.0, 0.3, 0.0, 0.0),
#'               c(0.0, 0.0, 0.1, 0.1))
#' 
#' matF <- rbind(c(0.0, 0.0, 5.0, 10.0),
#'               c(0.0, 0.0, 0.0,  0.0),
#'               c(0.0, 0.0, 0.0,  0.0),
#'               c(0.0, 0.0, 0.0,  0.0))
#' 
#' matC <- rbind(c(0.0, 0.0, 0.0, 0.0),
#'               c(0.0, 0.0, 0.0, 1.0),
#'               c(0.0, 2.0, 0.0, 0.0),
#'               c(0.0, 0.0, 0.0, 0.0))
#' 
#' dEntropy(matU, matC = matC)
#' dEntropy(matU, matF, matC, nSteps = 5)
#' dEntropy(matU, matF, matC, nSteps = 10)
#' dEntropy(matU, matF, matC, nSteps = 100)
#' @export dEntropy
dEntropy <- function(matU, matF, matC, startLife = 1, nSteps = 100) {
  
  # Error checks
  if (dim(matU)[1] != dim(matU)[2]) {
    stop("Matrix population model is not a square matrix")
  }
  if (!missing(matF) & !missing(matC)) {
    if (any(is.na(matF)) & any(is.na(matC))) {
      stop("NAs exist in both matF and matC")
    }
  }
  if (any(colSums(matU) > 1)) {
    warning("matU has at least one stage-specific survival value > 1")
  }

  # Create empty output list
  out <- NULL

  # Age-specific survivorship (lx)
  lx <- ageSpecificSurv(matU, startLife, nSteps)

  # Age-specific fertility (mx)
  if (!missing(matF)) {
    if (sum(matF, na.rm = TRUE) == 0) {
      warning("matF contains only zeros")
    }
    
    mx <- ageSpecificRepro(matU, matF, startLife, nSteps)
    out$Fec <- dEntropyCalc(lx, mx)
  }

  # Age-specific clonality (cx)
  if (!missing(matC)) {
    if (sum(matC, na.rm = TRUE) == 0) {
      warning("matC contains only zeros")
    }

    cx <- ageSpecificRepro(matU, matC, startLife, nSteps)
    out$Clo <- dEntropyCalc(lx, cx)
  }

  # Age-specific total reproduction (mxcx)
  if (!missing(matF) & !missing(matC)) {
    mxcx <- mx + cx
    out$FecClo <- dEntropyCalc(lx, mxcx)
  }
  
  return(out)
}


# Calculate Demetrius' entropy given lx and fx
dEntropyCalc <- function(lx, fx) {
  lxfx <- lx * fx

  # if lxfx == 0, log(lxfx) == -Inf; for entropy calc below, these -Inf can be
  #  converted to 0, because lim(x * log(x)) as x->0 is 0
  log_lxfx <- log(lxfx)
  log_lxfx[lxfx == 0] <- 0

  H <- -sum(lxfx * log_lxfx) / sum(lxfx)
  return(H)
}
