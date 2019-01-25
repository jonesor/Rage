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
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life.
#' @param N Maximum age to which age-specific traits will be calculated.
#' 
#' @return A vector of length \code{N+1}
#'   
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @seealso \code{\link{surv_conversion}}
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
#' mpm_to_lx(matU, startLife = 1, N = 20)
#' 
#' # age-specific survival probability
#' mpm_to_px(matU, startLife = 1, N = 20)
#' 
#' # age-specific mortality hazard
#' mpm_to_hx(matU, startLife = 1, N = 20)
#' 
#' # age-specific fecundity
#' mpm_to_mx(matU, matF, startLife = 1, N = 20)
#' 
#' @name age_from_stage
NULL


#' @rdname age_from_stage
#' @export mpm_to_lx
mpm_to_lx <- function(matU, startLife, N) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(startLife, matU)
  if (missing(N)) { stop("Argument N must be specified", call. = FALSE) }
  
  matUtemp <- matU
  lx <- vector(mode = 'numeric', length = N)
  
  for (i in 1:N) {
    lx[i] <- sum(matUtemp[,startLife])
    matUtemp <- matUtemp %*% matU
  }
  
  return(c(1, lx))
}


#' @rdname age_from_stage
#' @export mpm_to_px
mpm_to_px <- function(matU, startLife, N) {
  # leave argument validation to mpm_to_lx
  lx <- mpm_to_lx(matU, startLife, N)
  return(lx_to_px(lx))
}


#' @rdname age_from_stage
#' @export mpm_to_hx
mpm_to_hx <- function(matU, startLife, N) {
  # leave argument validation to mpm_to_lx
  lx <- mpm_to_lx(matU, startLife, N)
  return(lx_to_hx(lx))
}



#' @rdname age_from_stage
#' @export mpm_to_mx
mpm_to_mx <- function(matU, matR, startLife, N) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidMat(matR)
  checkValidStartLife(startLife, matU)
  if (missing(N)) { stop("Argument N must be specified", call. = FALSE) }
  
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
  return(c(0, mx))
}

