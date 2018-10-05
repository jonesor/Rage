#' Calculate age-specific survivorship from matrix population model
#'
#' This function uses age-from-stage decomposition methods to calculate
#' age-specific survivorship from a stage-classified matrix population model. A
#' detailed description of these methods can be found in section 5.3.1 of
#' Caswell (2001).
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life.
#' @param N Maximum age to which age-specific survivorship will be calculated.
#' @return A vector containing age-specific survivorship, from ages 0 through
#'   \code{N}.
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
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
#' ageSpecificSurv(matU, startLife = 1, N = 20)
#' 
#' @export ageSpecificSurv
ageSpecificSurv <- function(matU, startLife, N) {
  
  if (missing(startLife)) {
    stop("Argument startLife must be specified")
  }
  if (missing(N)) {
    stop("Argument N must be specified")
  }
  
  matUtemp <- matU
  lx <- vector(mode = 'numeric', length = N)
  
  for (i in 1:N) {
    lx[i] <- sum(matUtemp[,startLife])
    matUtemp <- matUtemp %*% matU
  }
  
  lx <- c(1, lx)
  return(lx)
}



#' Calculate age-specific reproduction from matrix population model
#'
#' This function uses age-from-stage decomposition methods to calculate
#' age-specific rates of reproduction from a stage-classified matrix population
#' model. A detailed description of these methods can be found in section 5.3.2
#' of Caswell (2001).
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life.
#' @param N Maximum age to which age-specific reproduction will be calculated.
#' @return A vector containing age-specific rates of reproduction, from ages 0
#'   through \code{N}.
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
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
#' ageSpecificRepro(matU, matF, startLife = 1, N = 20)
#' 
#' @export ageSpecificRepro
ageSpecificRepro <- function(matU, matR, startLife, N) {
  
  if (missing(startLife)) {
    stop("Argument startLife must be specified")
  }
  if (missing(N)) {
    stop("Argument N must be specified")
  }
  
  matUtemp <- matU
  mx <- vector(mode = 'numeric', length = N)
  
  for (i in 1:N) {
    # stageDist equivalent to: matUtemp %*% solve(diag(colSums(matUtemp)))
    stageDist <- apply(matUtemp, 2, function(x) x / sum(x))
    phi <- matR %*% stageDist
    mx[i] <- sum(phi[,startLife])
    matUtemp <- matUtemp %*% matU
  }
  
  mx[is.nan(mx)] <- 0
  mx <- c(0, mx)
  return(mx)
}

