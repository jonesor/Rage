#' Rearrange stages of a matrix population model to segregate reproductive and
#' non-reproductive stages
#' 
#' Rearrange stages of a matrix population model so that all inter-reproductive
#' stages fall in the final rows/columns of the matrix. This is a preparatory
#' step to collapsing the matrix model into a standardized set of stages (e.g.
#' propagule, pre-reproductive, reproductive, and post-reproductive).
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matF The sexual component of a matrix population model (i.e. a square
#'   projection matrix reflecting transitions due to sexual reproduction)
#' @param matC The clonal component of a matrix population model (i.e. a square
#'   projection matrix reflecting transitions due to clonal reproduction).
#'   Defaults to \code{NULL}, indicating no clonal reproduction (i.e.
#'   \code{matC} is a matrix of zeros).
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
#' matU <- rbind(c(0.1,   0,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0,   0),
#'               c(  0, 0.3, 0.3, 0.1,   0),
#'               c(  0,   0, 0.4, 0.4, 0.1),
#'               c(  0,   0,   0, 0.1, 0.4))
#' 
#' matF <- rbind(c(  0, 1.1,   0, 1.6,   0),
#'               c(  0, 0.8,   0, 0.4,   0),
#'               c(  0,   0,   0,   0,   0),
#'               c(  0,   0,   0,   0,   0),
#'               c(  0,   0,   0,   0,   0))
#' 
#' reproStages <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
#' matrixStages <- c('prop', 'active', 'active', 'active', 'active')
#' 
#' rearrangeMatrix(matU, matF, reproStages = reproStages,
#'                 matrixStages = matrixStages)
#' 
#' @export rearrangeMatrix
rearrangeMatrix <- function(matU, matF, matC = NULL, reproStages,
                            matrixStages) {

  # validate arguments
  checkValidMat(matU)
  checkValidMat(matF)
  if (!is.null(matC)) checkValidMat(matC, warn_all_zero = FALSE)
  if (ncol(matU) != ncol(matF) ||
        ncol(matU) != length(reproStages) ||
          length(reproStages) != length(matrixStages)) {
    stop("Arguments do not correspond to MPM of single dimension",
         call. = FALSE)
  }
  
  # populate matC with zeros, if NULL
  if (is.null(matC)) {
    matC <- matrix(0, nrow = ncol(matF), ncol = ncol(matF))
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
