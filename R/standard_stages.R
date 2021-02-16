#' Group the stages of a matrix population model into a standardized set of
#' stage classes
#'
#'@description 
#' Group the stages of a matrix population model into a standardized set of
#' stage classes. Matrix population model have stages representing different
#' moments of the life cycle. Standard stages allows the user to obtain a vector 
#' describing the position of the different stages in the matrix population
#' model; e.g. to identify which are the reproductive stages in a matrix 
#' population model. This is particularly useful to use in combination with 
#' functions that require entering some particular stages, see
#' \code{\link{net_repo_rate}}. 
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
#' @param includeProp (logical) include propagule stage. default: \code{TRUE}.
#' if \code{TRUE}, propagule stage (if present) is given back in result. If
#' \code{FALSE}, it's included into the pre-reproductive stage
#' @param includePost (logical) include post-reproductive stage. default:
#' \code{TRUE}. if \code{TRUE}, post-reproductive stage (if present) is given
#' back in result. If \code{FALSE}, it's included into the reproductive
#' stage
#' 
#' @return A list with four potential elements:
#'   \item{propStages}{Position of the propagule stages}
#'   \item{preRepStages}{Position of the reproductive stages}
#'   \item{matF}{Sexual reproduction component of the collapsed projection
#'   matrix}
#'   \item{matC}{Clonal reproduction component of the collapsed projection
#'   matrix}
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
#' # combine post-reproductive and reproductive
#' standard_stages(matF, reproStages, matrixStages, includePost = FALSE)
#' 
#' @export standard_stages
standard_stages <- function(matF, reproStages, matrixStages) {
  # validate arguments
  checkValidMat(matF, warn_all_zero = FALSE)
  if (ncol(matF) != length(reproStages) ||
      length(reproStages) != length(matrixStages)) {
    stop("Arguments do not correspond to MPM of single dimension",
         call. = FALSE)
  }
  if (!any(reproStages == TRUE)) {
    stop("Cannot identify standardized stages because no stages are ",
         "reproductive (i.e. at least one element of reproStages must be TRUE)",
         call. = FALSE)
  }
  
  
  if ("prop" %in% matrixStages) {
    propStage <- which(matrixStages == "prop")
    # if (!includeProp) {
    #   warning("Propagule stage exists, but includeProp is set to FALSE")
    # }
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
