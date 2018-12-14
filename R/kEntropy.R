#' Calculate Keyfitz's entropy from a matrix population model
#' 
#' This function calculates Keyfitz's entropy from a matrix population model, by
#' first using age-from-stage decomposition methods to estimate age-specific
#' survivorship (lx).
#' 
<<<<<<< HEAD
#' @param matU A matrix containing only survival-dependent processes (e.g. progression,
#' stasis, retrogression).
#' @param startLife The first stage at which the author considers the beginning
#' of life in the life cycle of the species. It defaults to the first stage.
#' @param nSteps A cutoff for the decomposition of age-specific survival ('lx'), and when pertinent,
#' for age-specific sexual reproduction ('mx') and age-specific clonal reproduction ('cx'). This allows
#' excluding mortality and fertility plateaus. See function 'qsdConverge' for more information. When not
#' specified, this argument defaults to 100.
#' @param trapeze A logical argument indicating whether the trapezoidal
#' approximation should be used for approximating the definite integral.
#' @param trimlx A numerical argument determining how lx should be trimmed. 
#' For example, a value of 0.05 means that lx values are retained only if they 
#' are greater than or equal to 0.05. NB: 0 values in lx mean that Keyfitz 
#' entropy cannot be estimated.
#' @return Returns an estimate of Keyfitz' life table entropy based on an lx
#' (survivorship) vector obtained from matU
=======
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param nSteps The age-cutoff for the decomposition of age-specific survival
#'   (lx). This allows the user to exclude ages after which mortality or
#'   fertility has plateaued (see function \code{\link{qsdConverge}} for more
#'   information). Defaults to 100.
#' @param trapeze A logical argument indicating whether the composite trapezoid
#'   approximation should be used for approximating the definite integral.
#' @return Returns an estimate of Keyfitz's life table entropy.
>>>>>>> devel
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @references Keyfitz, N. (1977) Applied Mathematical Demography. New York:
#'   Wiley.
#' @examples
#' matU <- rbind(c(0.2,   0,   0,   0),
#'               c(0.3, 0.4, 0.1,   0),
#'               c(0.1, 0.1, 0.2, 0.3),
#'               c(  0, 0.2, 0.6, 0.5))
#' 
#' kEntropy(matU, nSteps = 10)
#' kEntropy(matU, nSteps = 20)
#' kEntropy(matU, nSteps = 100)
#' kEntropy(matU, nSteps = 100, trapeze = TRUE)
#' 
#' @export kEntropy
<<<<<<< HEAD
kEntropy <- function(matU, startLife = 1, nSteps = 1000, trapeze = FALSE, trimlx = 0.05){
=======
kEntropy <- function(matU, startLife = 1, nSteps = 100, trapeze = FALSE) {
>>>>>>> devel
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(startLife, matU)
  
  # Age-specific survivorship (lx)
  lx <- ageSpecificSurv(matU, startLife, nSteps)
  lx <- lx[1:max(which(lx > 0))] # remove ages at/beyond which lx is 0 or NA
  
<<<<<<< HEAD
  #Trim lx
  lx <- lx[lx>=trimlx]
  
  if(trapeze == TRUE){
    ma <- function(x,n=2){stats::filter(x,rep(1/n,n), sides=2)}
    lx2 <- na.omit(as.vector(ma(lx)))
    return(-sum(lx2*log(lx2))/sum(lx2))
    }else{
      return(-sum(lx*log(lx))/sum(lx))
    }
=======
  # Calculate Keyfitz's entropy
  if (trapeze == TRUE) {
    H <- -TrapezoidRule(lx * log(lx)) / TrapezoidRule(lx)
  } else {
    H <- -sum(lx * log(lx)) / sum(lx)
>>>>>>> devel
  }
  
  return(H) 
}


# Composite Trapezoid Rule for approximating a definite integral;
# modified to work with discrete y values that start at x = 0, and have
# h = delta_x = 1 (i.e. lx)
TrapezoidRule <- function(y) {
  n <- length(y)
  h <- 1   # h = (b - a) / n  = (length(y) - 0) / n
  return((h / 2) * (y[1] + 2 * sum(y[2:(n-1)]) + y[n]))
}
