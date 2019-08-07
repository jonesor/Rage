#' Calculate age-specific traits from a from matrix population model
#'
#' These functions use age-from-stage decomposition methods to calculate
#' age-specific survivorship (lx), survival probability (px), mortality hazard
#' (hx), or reproduction (mx) from a matrix population model. A detailed
#' description of these methods can be found in sections 5.3.1 and 5.3.2 of
#' Caswell (2001).
#' 
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @param start The index of the first stage at which the author considers the
#'   beginning of life. Defaults to \code{1L}.
#' @param xmax Maximum age to which age-specific traits will be calculated
#'   (defaults to \code{1e5}).
#' @param lxCrit Minimum value of lx to which age-specific traits will be
#'   calculated (defaults to \code{1e-4}).
#' 
#' @return A vector
#' 
#' @note The output vector is calculated recursively until the age class (x)
#'   reaches \code{xmax} or survivorship (lx) falls below \code{lxCrit} —
#'   whichever comes first. To force calculation to \code{xmax}, set
#'   \code{lxCrit = 0}. Conversely, to force calculation to \code{lxCrit}, set
#'   \code{xmax = Inf}.
#'   
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @seealso \code{\link{lifetable_convert}}
#' 
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' 
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
#' # age-specific survivorship
#' mpm_to_lx(matU, start = 1, xmax = 20)
#' 
#' # age-specific survival probability
#' mpm_to_px(matU, start = 1, xmax = 20)
#' 
#' # age-specific mortality hazard
#' mpm_to_hx(matU, start = 1, xmax = 20)
#' 
#' # age-specific fecundity
#' mpm_to_mx(matU, matF, start = 1, xmax = 20)
#' 
#' @name age_from_stage
NULL


matU <- rbind(c(0.0, 0.0, 0.0, 0.0),
              c(0.3, 0.3, 0.1, 0.0),
              c(0.0, 0.2, 0.1, 0.6),
              c(0.0, 0.1, 0.1, 0.3))

matR <- rbind(c(0.0, 0.0, 1.4, 2.1),
              c(0.0, 0.0, 0.0, 0.0),
              c(0.0, 0.0, 0.0, 0.0),
              c(0.0, 0.0, 0.0, 0.0))


#' @rdname age_from_stage
#' @export mpm_to_lx
mpm_to_lx <- function(matU, start = 1L, xmax = 1e5, lxCrit = 1e-4) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU)
  
  lx <- 1.0
  lx_vec <- lx
  n <- rep(0.0, nrow(matU))
  n[start] <- 1.0
  t <- 0L
  
  while (lx > lxCrit & t < xmax) {
    n <- matU %*% n
    lx <- sum(n)
    t <- t + 1L
    lx_vec[t] <- lx
  }
  
  return(c(1.0, lx_vec))
}


#' @rdname age_from_stage
#' @export mpm_to_px
mpm_to_px <- function(matU, start = 1L, xmax = 1e5, lxCrit = 1e-4) {
  # leave argument validation to mpm_to_lx
  lx <- mpm_to_lx(matU, start, xmax, lxCrit)
  return(lx_to_px(lx))
}


#' @rdname age_from_stage
#' @export mpm_to_hx
mpm_to_hx <- function(matU, start = 1L, xmax = 1e5, lxCrit = 1e-4) {
  # leave argument validation to mpm_to_lx
  lx <- mpm_to_lx(matU, start, xmax, lxCrit)
  return(lx_to_hx(lx))
}



#' @rdname age_from_stage
#' @export mpm_to_mx
mpm_to_mx <- function(matU, matR, start = 1L, xmax = 1e5, lxCrit = 1e-4) {
  
  # validate arguments (leave rest to mpm_to_lx)
  checkValidMat(matR)
  
  N <- length(mpm_to_lx(matU, start, xmax, lxCrit))
  
  tempU <- matU
  mx <- vector(mode = 'numeric', length = N)
  
  for (i in 1:N) {
    # stageDist equivalent to: matUtemp %*% solve(diag(colSums(matUtemp)))
    stageDist <- apply(tempU, 2, function(x) x / sum(x))
    phi <- matR %*% stageDist
    mx[i] <- sum(phi[,start])
    tempU <- tempU %*% matU
  }
  
  mx[is.nan(mx)] <- 0
  return(mx)
}
