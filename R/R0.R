#' Calculate net reproductive value from a matrix population model
#'
#' Calculate net reproductive value (i.e. the per-generation population growth
#' rate) from a matrix population model
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @return Returns the net reproductive value.  If \code{matU} is singular
#'   (often indicating infinite life expectancy), returns \code{NA}.
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
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
#' R0(matU, matF)
#' 
#' @importFrom popbio lambda
#' @export R0
R0 <- function(matU, matR) {
  
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
                 
  # matrix dimensions
  matDim <- nrow(matU)
  
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)
  
  # calculate R0
  # first check for singular matU (if singular, R0 = NA)
  if (class(N) == 'try-error' && grepl('singular', N[1])) {
    R0 <- NA_real_
  } else {
    R <- matR %*% N
    R0 <- popbio::lambda(R)
  }
  
  return(R0)
}
