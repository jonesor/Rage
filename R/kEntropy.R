#' Calculate Keyfitz's entropy
#' 
#' This function calculates Keyfitz's entropy from a matrix population model, by
#' first using age-from-stage decomposition methods to estimate age-specific
#' survivorship (lx).
#' 
#' @param matU A square matrix containing only survival-related transitions
#'   (i.e. progression, stasis, retrogression).
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param nSteps The age-cutoff for the decomposition of age-specific survival
#'   (lx). This allows the user to exclude ages after which mortality or
#'   fertility has plateaued (see function \code{qsdConverge} for more
#'   information). Defaults to 100.
#' @param trapeze A logical argument indicating whether the trapezoidal
#'   approximation should be used for approximating the definite integral.
#' @return Returns an estimate of Keyfitz's life table entropy.
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @references Keyfitz, N. (1977) Applied Mathematical Demography. New York:
#'   Wiley.
#' @examples
#' matU <- rbind(c(0.2, 0.0, 0.0, 0.0),
#'               c(0.3, 0.4, 0.1, 0.0),
#'               c(0.1, 0.1, 0.2, 0.3),
#'               c(0.0, 0.2, 0.6, 0.5))
#' 
#' kEntropy(matU, nSteps = 10)
#' kEntropy(matU, nSteps = 20)
#' kEntropy(matU, nSteps = 100)
#' kEntropy(matU, nSteps = 100, trapeze = TRUE)
#' @export kEntropy
kEntropy <- function(matU, startLife = 1, nSteps = 100, trapeze = FALSE) {
  
  # Error checks
  if (dim(matU)[1] != dim(matU)[2]) {
    stop("Matrix population model is not a square matrix")
  } 
  if (any(is.na(matU))) {
    stop("NAs exist in matU")
  } 
  if (any(colSums(matU) > 1)) {
    warning("matU has at least one stage-specific survival value > 1")
  } 
  
  # Age-specific survivorship (lx)
  lx <- ageSpecificSurv(matU, startLife, nSteps)
  lx <- lx[1:max(which(lx > 0))] # remove ages at/beyond which lx is 0 or NA
  
  # Calculate Keyfitz's entropy
  if (trapeze == TRUE) {
    H <- -Trapezoid(lx * log(lx)) / Trapezoid(lx)
  } else {
    H <- -sum(lx * log(lx)) / sum(lx)
  }
  
  return(H) 
}


Trapezoid <- function(y) {
  n <- b <- length(y)
  a <- 0
  h <- (b - a) / n
  
  out <- (h / 2) * (y[1] + 2 * sum(y[2:(n-1)]) + y[n])
  return(out)
}
