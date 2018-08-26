#' Calculate net reproductive value
#'
#' Calculate net reproductive value from a matrix population model.
#'
#' @param matU A matrix containing only survival-related transitions (i.e.
#'   progression, stasis, retrogression).
#' @param matR A matrix containing only reproduction-related transitions (either
#'   sexual, clonal, or both; i.e. \code{matF}, \code{matC}, or \code{mat +
#'   matC}).
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @return Returns the net reproductive value.
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' @examples
#' matU <- rbind(c(0.0, 0.0, 0.0, 0.0),
#'               c(0.5, 0.0, 0.0, 0.0),
#'               c(0.0, 0.3, 0.0, 0.0),
#'               c(0.0, 0.0, 0.1, 0.1))
#' 
#' matF <- rbind(c(0, 0, 5, 9),
#'               c(0, 0, 0, 0),
#'               c(0, 0, 0, 0),
#'               c(0, 0, 0, 0))
#' 
#' R0(matU, matF)
#' R0(matU, matF, startLife = 2)
#' @export R0
R0 <- function(matU, matR, startLife = 1) {
  
  # Validate inputs
  if (nrow(matU) != ncol(matU)) {
    stop("matU is not a square matrix")
  }
  if (nrow(matR) != ncol(matR)) {
    stop("matU is not a square matrix")
  }
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
  
  # Fundamental matrix, which states the amount of time units spent in each
  #  stage on average
  N <- solve(diag(matDim) - matU)
  
  # calculate R0
  R0_matR <- matR %*% N
  R0 <- R0_matR[startLife, startLife]
  
  return(R0)
}
