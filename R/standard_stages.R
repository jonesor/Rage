#' Identify stages corresponding to different parts of the reproductive life
#' cycle
#'
#'@description 
#' Identify the stages of a matrix population model that correspond to
#' different parts of the reproductive life cycle, namely propagule,
#' pre-reproductive, reproductive and post-reproductive. These classifications
#' are used to standardise matrices to allow comparisons across species with
#' different life cycle structures, see \code{\link{mpm_standardize}}.
#' 
#' Assumes that fecundity and mean fecundity matrices have been 
#' rearranged so that non-reproductive stages are in the final rows/columns.
#' Output indicates groupings to be used when collapsing the matrix model.
#'
#' @param matF The sexual component of a matrix population model (i.e. a square
#'   projection matrix reflecting transitions only due to \emph{sexual} 
#'   reproduction). It assumes that it has been rearranged so that 
#'   non-reproductive stages are in the final rows/columns.
#' @param reproStages Logical vector identifying which stages are reproductive.
#' @param matrixStages (character) vector of stages, values are "prop"
#' (propagule), "active", and "dorm" (dormant).
#' 
#' @return A list with four elements:
#'   \item{propStages}{Position of the propagule stages}
#'   \item{preRepStages}{Position of the pre-reproductive stages}
#'   \item{repStages}{Position of the reproductive stages}
#'   \item{postRepStages}{Position of the post-reproductive stages}
#' 
#' @author Rob Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @note Dormant stages are not currently handled.
#' @seealso \code{\link{mpm_standardize}}
#' 
#' @examples
#' matF <- rbind(c(  0, 1.1,   0, 1.6,   0),
#'               c(  0, 0.8,   0, 0.4,   0),
#'               c(  0,   0,   0,   0,   0),
#'               c(  0,   0,   0,   0,   0),
#'               c(  0,   0,   0,   0,   0))
#'
#' reproStages <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
#' matrixStages <- c('prop', 'active', 'active', 'active', 'active')
#' standard_stages(matF, reproStages, matrixStages)
#' 
#' @export standard_stages
standard_stages <- function(matF, reproStages, matrixStages) {

  # validate arguments
  checkValidMat(matF, warn_all_zero = FALSE)
  if (ncol(matF) != length(reproStages) ||
      length(reproStages) != length(matrixStages)) {
    stop("Arguments do not correspond to MPM of the same dimension",
         call. = FALSE)
  }
  if (!any(reproStages == TRUE)) {
    stop("Cannot identify standardised stages because no stages are ",
         "reproductive (i.e. at least one element of reproStages must be TRUE)",
         call. = FALSE)
  }
  
  if ("prop" %in% matrixStages) {
    propStage <- which(matrixStages == "prop")
  } else {
    propStage <- NA
  }
  
  # set max reproductive stage
  maxRep <- max(which(reproStages == TRUE))
  
  # prerep
  matDim <- nrow(matF)
  Rep <- which(reproStages == TRUE)
  
  if (min(Rep) == 1) {
    preRep <- NA
  } else if (!is.na(propStage[1]) && (min(Rep) - max(propStage) == 1)) {
    preRep <- NA
  } else {
    preRep <- min(which(matrixStages == "active")):(min(Rep) - 1)
  }
  
  # postrep
  if (maxRep == matDim) {
    postRep <- NA
  } else {
    postRep <- (maxRep + 1):matDim
  }
  
  return(list(propStages = propStage,
              preRepStages = preRep,
              repStages = Rep,
              postRepStages = postRep))
}
