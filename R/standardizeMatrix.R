#' Transform a matrix population model to a standardized form
#'
#' Transform a matrix population model to a standardized set of stage classes
#' (e.g. propagule, pre-reproductive, reproductive, and post-reproductive). The
#' transition rates in the standardized matrix are a weighted average of the
#' transition rates from the relevant stages of the original matrix, weighted by
#' the relative proportion of each stage class expected at the stable
#' distribution.
#' 
#' @param matU A square matrix containing only survival-related transitions
#'   (i.e. progression, stasis, retrogression).
#' @param matF A square matrix containing only sexual reproduction-related
#'   transitions.
#' @param matC A square matrix containing only clonal reproduction-related
#'   transitions. The default is \code{NULL}, indicating no clonal
#'   reproduction (i.e. \code{matC} is a matrix of zeros).
#' @param reproStages Logical vector indicating which stages are reproductive
#' @param matrixStages Character vector of matrix stage types (e.g. "propagule",
#'   "active", or "dormant")
#' @return A list with four elements reflecting the standardized matrix, and
#'   it's components:
#'   \item{matA}{Standardized projection matrix}
#'   \item{matU}{Survival component of the standardized projection matrix}
#'   \item{matF}{Sexual reproduction component of the standardized projection matrix}
#'   \item{matC}{Clonal reproduction component of the standardized projection matrix}
#' @details This function is a wrapper for the functions
#'   \code{\link{rearrangeMatrix}}, \code{\link{reprodStages}} and
#'   \code{\link{collapseMatrix}}, which it calls in sequence.
#' @note The method used by this function to collapse a matrix population model
#'   preserves the equilibrium population growth rate (\eqn{lamda}) and relative
#'   stable distribution, but is not expected to preserve other traits such as
#'   relative reproductive values, sensitivities, net reproductive rates, life
#'   expectancy, etc.
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @examples
#' matU <- rbind(c(  0,   0,    0,    0),
#'               c(0.5,   0,    0,    0),
#'               c(  0, 0.3,    0,    0),
#'               c(  0,   0,  0.2,  0.1))
#' 
#' matF <- rbind(c(  0,   0,  1.1,  1.6),
#'               c(  0,   0,  0.8,  0.4),
#'               c(  0,   0,    0,    0),
#'               c(  0,   0,    0,    0))
#'
#' matC <- rbind(c(  0,   0,  0.4,  0.5),
#'               c(  0,   0,  0.3,  0.1),
#'               c(  0,   0,    0,    0),
#'               c(  0,   0,    0,    0))
#'
#' reproStages <- c(FALSE, FALSE, TRUE, TRUE)
#' matrixStages <- c('prop', 'active', 'active', 'active')
#'
#' standardizeMatrix(matU, matF, matC, reproStages, matrixStages)
#' 
#' @export standardizeMatrix
standardizeMatrix <- function(matU, matF, matC = NULL, reproStages,
                              matrixStages) {
  
  # populate matC NULL, populate with zeros
  if (is.null(matC)) {
    matC <- matrix(0, nrow = ncol(matF), ncol = ncol(matF))
  }
  
  # put non-reproductive stages at the end of the matrix
  rearr <- rearrangeMatrix(matU, matF, matC, reproStages, matrixStages)
  
  # defines which columns need to be collapsed for each of the four stages
  collapse <- reprodStages(rearr$matF,
                           rearr$reproStages,
                           rearr$matrixStages)
  
  # collapse
  out <- collapseMatrix(matU, matF, matC, collapse = collapse)
  
  # NA out of the standardized stages that were non-existant
  out$matA[is.na(collapse), is.na(collapse)] <- NA_real_
  out$matU[is.na(collapse), is.na(collapse)] <- NA_real_
  out$matF[is.na(collapse), is.na(collapse)] <- NA_real_
  out$matC[is.na(collapse), is.na(collapse)] <- NA_real_

  return(out)
}
