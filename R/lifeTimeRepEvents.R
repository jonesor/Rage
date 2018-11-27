#' Calculate the timing of lifetime reproductive events
#'
#' Applies Markov chain approaches to decompose various moments along the
#' age-specific trajectory of reproduction of individuals in a matrix population
#' model, including the probability of achieving maturity, age at first
#' reproduction, mean life expectancy conditional on maturity, and life
#' expectancy for mature individuals.
#' 
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @return A list with four elements:
#' \item{p}{probability of achieving reproductive maturity}
#' \item{La}{mean age at maturity (in the same units as the matrix population
#' model sampling periodicity)}
#' \item{meanLifeExpectancy}{mean life expectancy conditional on entering the
#' life cycle at the value of \code{startLife}}
#' \item{remainingMatureLifeExpectancy}{Life expectancy from mean maturity,
#' calculated as \code{meanLifeExpectancy - La}. This value can be negative
#' because both mean life expectancy and mean age at maturity are means of their
#' respective distributions, and their distributions can indeed overlap}
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <hcaswell@@whoi.edu>
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
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
#' lifeTimeRepEvents(matU, matF, startLife = 1)
#' 
#' @importFrom MASS ginv
#' @export lifeTimeRepEvents
lifeTimeRepEvents <- function(matU, matR, startLife = 1) {

  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidMat(matR)
  checkValidStartLife(startLife, matU)
  
  # initialize output list and get matrix dimension
  out <- list()
  matDim <- nrow(matU)
  
  # Which stages reproductive?
  fecLifeStages <- apply(matR, 2, function(x) sum(x, na.rm = TRUE) > 0)
  
  # Probability of survival to first sexual reproductive event
  # Note: U matrices are called 'T' in Caswell (2001)
  Uprime <- matU
  Uprime[,fecLifeStages] <- 0
  
  Mprime <- matrix(0, nrow = 2, ncol = matDim)
  
  for (p in 1:matDim) {
    if (fecLifeStages[p]) {
      Mprime[2, p] <- 1
    } else {
      Mprime[1, p] <- 1 - colSums(matU)[p]
    }
  }
  
  Bprime <- Mprime %*% (MASS::ginv(diag(matDim) - Uprime))
  out$p <- Bprime[2, startLife]
  
  # Age at first reproduction (La; Caswell 2001, p 124)
  D <- diag(c(Bprime[2,]))
  Uprimecond <- D %*% Uprime %*% MASS::ginv(D)
  expTimeReprod <- colSums(MASS::ginv(diag(matDim) - Uprimecond))
  out$La <- expTimeReprod[startLife]
  
  # Life expectancy conditional on entering the life cycle in the first
  #  reproductive stage
  firstFecLifeStage <- min(which(fecLifeStages))
  N <- solve(diag(matDim) - matU)
  out$meanLifeExpectancy <- colSums(N)[firstFecLifeStage]
  
  # Life expectancy from mean age of first reproduction maturity
  out$remainingMatureLifeExpectancy <- colSums(N)[startLife] - out$La
  
  return(out)
}
