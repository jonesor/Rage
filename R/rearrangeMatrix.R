#' Rearrange stages of a matrix population model to segregate reproductive and
#' non-reproductive stages
#' 
#' Rearrange stages of a matrix population model so that all inter-reproductive
#' stages fall in the final rows/columns of the matrix. This is a preparatory
#' step to collapsing the matrix model into a standardized set of stages (e.g.
#' propagule, pre-reproductive, reproductive, and post-reproductive).
#'
#' @param matU A square matrix containing only survival-related transitions
#'   (i.e. progression, stasis, retrogression).
#' @param matF A square matrix containing only sexual reproduction-related
#'   transitions.
#' @param matC A square matrix containing only clonal reproduction-related
#'   transitions.
#' @param matrixStages A character vector identifying organized matrix stages
#' @param reproStages Logical vector identifying which stages reproductive
#' @return Returns a list with 6 elements:
#' \item{matU}{Rearranged survival matrix}
#' \item{matF}{Rearranged sexual reproduction matrix}
#' \item{matC}{Rearranged clonal reproduction matrix}
#' \item{matrixStages}{Rearranged vector of organized matrix stages}
#' \item{reproStages}{Rearranged logical vector of reproductive stages}
#' \item{nonRepInterRep}{Numeric index for any rearranged inter-reproductive
#'  stages}
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @seealso \code{\link{standardizeMatrix}}
#' @examples
#' matU <- rbind(c(0, 0, 0, 0, 0), c(0.1, 0.16, 0, 0, 0), c(0.2, 0.23, 0.12, 0,
#' 0), c(0, 0, 0.34, 0.53, 0), c(0, 0, 0, 0.34, 0))
#' 
#' matF <- rbind(c(0, 0, 0, 0, 0), c(0, 0.2, 0, 0.1, 0), c(0, 0.2, 0, 0.1, 0),
#' c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
#' 
#' reproStages <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
#' matrixStages <- c('prop', 'active', 'active', 'active', 'active')
#' rearrangeMatrix(matU, matF, reproStages = reproStages,
#'                 matrixStages = matrixStages)
#' @export rearrangeMatrix
rearrangeMatrix <- function(matU, matF, matC = NULL, reproStages,
                            matrixStages) {

  # populate matC with zeros, if NULL
  if (is.null(matC)) {
    matC <- matrix(0, nrow = ncol(matF), ncol = ncol(matF))
  }
  
  # validate arguments
  if (ncol(matU) != ncol(matF) |
      ncol(matU) != length(reproStages) |
      length(reproStages) != length(matrixStages)) {
    stop("Arguments do not correspond to MPM of single dimension")
  }
  
  # preliminaries
  matDim <- ncol(matF)              # matrix dim
  Rep <- which(reproStages == TRUE) # reproductive stage indices
  out <- NULL                       # initalize output list
  
  # which stages reproductive, incl. inter-reproductive
  if (length(Rep) > 0) {
    allRep <- Rep[1]:Rep[length(Rep)]
  } else {
    allRep <- integer(0)
  }
  
  # which stages non-reproductive inter-reproductive
  nonRepInterRep <- allRep[which(!allRep %in% Rep)]
  
  if (length(nonRepInterRep) > 0) {
    # if there are non-repro inter-repro stages, rearrange
    allElseStages <- which(!(1:matDim %in% nonRepInterRep))
    rearrange <- c(allElseStages, nonRepInterRep)
    out$matU <- matU[rearrange, rearrange]
    out$matF <- matF[rearrange, rearrange]
    out$matC <- matC[rearrange, rearrange]
    out$matrixStages <- matrixStages[rearrange]
    out$reproStages <- reproStages[rearrange]
  } else {
    # else, no need to rearrange
    out$matU <- matU
    out$matF <- matF
    out$matC <- matC
    out$matrixStages <- matrixStages
    out$reproStages <- reproStages
  }
  
  # inter-repro stages that were moved to the last column(s)
  out$nonRepInterRep <- ifelse(length(nonRepInterRep) > 0, nonRepInterRep, NA)
  
  return(out)
}
