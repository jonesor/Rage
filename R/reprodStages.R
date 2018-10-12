#' Group the stages of a matrix population model into a standardized set of
#' stage classes
#'
#' Assumes that fecundity and mean fecundity matrices have been rearranged so
#' that non-reproductive stages are in the final rows/columns.
#'
#' Output indicates groupings to be used when collapsing the matrix model.
#'
#' @param matF The sexual component of a matrix population model (i.e. a square
#'   projection matrix reflecting transitions due to sexual reproduction) —
#'   rearranged so that non-reproductive stages are in the final rows/columns
#' @param reproStages Logical vector identifying which stages reproductive
#' @param matrixStages (character) vector of stages, values are "prop"
#' (propagule), "active", and "dorm" (dormant)
#' @param includeProp (logical) include propagule stage. default: \code{TRUE}.
#' if \code{TRUE}, propagule stage (if present) is given back in result. If
#' \code{FALSE}, it's included into the pre-reproductive stage
#' @param includePost (logical) include post-reproductive stage. default:
#' \code{TRUE}. if \code{TRUE}, post-reproductive stage (if present) is given
#' back in result. If \code{FALSE}, it's included into the reproductive
#' stage
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @note Dormant stages are not currently handled.
#' @seealso \code{\link{standardizeMatrix}}
#' @examples
#' matF <- rbind(c(  0, 1.1,   0, 1.6,   0),
#'               c(  0, 0.8,   0, 0.4,   0),
#'               c(  0,   0,   0,   0,   0),
#'               c(  0,   0,   0,   0,   0),
#'               c(  0,   0,   0,   0,   0))
#'
#' reproStages <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
#' matrixStages <- c('prop', 'active', 'active', 'active', 'active')
#' reprodStages(matF, reproStages, matrixStages)
#' 
#' # combine post-reproductive and reproductive
#' reprodStages(matF, reproStages, matrixStages, includePost = FALSE)
#' @export reprodStages
reprodStages <- function(matF, reproStages, matrixStages, includeProp = TRUE,
                         includePost = TRUE) {
  
  # FIXME: can we combine propagule and pre-reproductive?
  # FIXME: what about dormant stages?
  # FIXME: once above solved, then fix logic internally
  # FIXME: Scott and Tamora added extra parameter flags for switching whether
  # to include propagule and post-reproductive stages but these are not used as
  # yet. On reflection this is probably going too far and we should let the
  # matrix structure guide the process. However, in animals we may wish to
  # combine rep and post-rep to obtain a two-stage matrix model.
  
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
