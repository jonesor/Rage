#' Group the stages of a matrix model into propagule, pre-reproductive,
#' reproductive and post-reproductive stages.
#'
#' Assumes that fecundity and mean fecundity matrices have been rearranged so
#' that non-reproductive stages are in the final rows/columns.
#'
#' Output indicates groupings to be used when collapsing the matrix model.
#'
#' Note: dormant stages not handled.
#'
#' FIXME: ROB says possibly put output as a named list instead of vector
#' FIXME: can we combine propagule and pre-reproductive?
#' FIXME: what about dormant stages?
#' FIXME: once above solved, then fix logic internally
#'
#' FIXME: Scott and Tamora added the extra parameter flags for switching whether
#' to include propagule and post-reproductive stages but these are not used as
#' yet. On reflection this is probably going too far and we should let the
#' matrix structure guide the process. However, in animals we may wish to
#' combine rep and post-rep to obtain a two-stage matrix model.
#'
#' @export
#' @param matF (matrix) fecundity matrix, rearranged so that non-reproductive
#' stages are in the final rows/columns
#' @param post (integer) indicator for post-reproductive stages
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
#' @examples
#' matU <- rbind( c(0, 0, 0, 0, 0), c(0.18, 0.16, 0, 0, 0), c(0.29, 0.23, 0.12,
#' 0, 0), c(0, 0, 0.34, 0.53, 0), c(0, 0, 0, 0.87, 0) )
#'
#' matF <- rbind( c(0, 0.13, 0.96, 0, 0), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0),
#' c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0) )
#'
#' reproStages <- c(FALSE, TRUE, TRUE, FALSE, FALSE)
#' matrixStages <- c('prop', 'active', 'active', 'active', 'active')
#' reprodStages(matF, post = 5, reproStages, matrixStages)
#' 
#' # combine post-reproductive and reproductive
#' reprodStages(matF, post = 5, reproStages, matrixStages, includePost = F)
#' \dontrun{
#' ## NOT ALLOWED
#' reprodStages(matF, post = 5, reproStages, matrixStages, c('prop', 'rep', 'postrep'))
#' reprodStages(matF, post = 5, reproStages, matrixStages, c('prop', 'postrep'))
#' reprodStages(matF, post = 5, reproStages, matrixStages, c('prerep', 'postrep'))
#' reprodStages(matF, post = 5, reproStages, matrixStages, c('prop', 'prerep'))
#' }
reprodStages <- function(matF, post, reproStages, matrixStages,
                         includeProp = TRUE, includePost = TRUE) {
  
  if ("prop" %in% matrixStages) {
    propStage <- which(matrixStages == "prop")
    if (!includeProp) {
      warning("Propagule stage exists, but includeProp is set to FALSE")
    }
  } else {
    propStage <- NA
  }
  
  # set max reproductive stage
  maxRep <- which.max(reproStages == 1)
  
  # prerep
  matDim <- dim(matF)[1]
  Rep <- which(reproStages == 1)
  if (length(Rep) > 0) {
    if (min(Rep) == 1) {
      preRep <- NA
    } else if (!is.na(propStage[1]) && (min(Rep) - max(propStage) == 1)) {
      preRep <- NA
    } else {
      preRep <- min(which(matrixStages == "active")):(min(Rep) - 1)
    }
  } else {
    preRep <- NA
  }
  
  # postrep
  if (length(post) == 0 && maxRep == matDim) {
    postRep <- NA
  } else {
    postRep <- (maxRep + 1):matDim
  }
  
  if (length(propStage) > 1) {
    propStages <- paste(propStage[1], "-", propStage[length(propStage)], sep = "")
  } else {
    propStages <- as.character(propStage)
  }
  if (length(preRep) > 1) {
    preRepStages <- paste(preRep[1], "-", preRep[length(preRep)], sep = "")
  } else {
    preRepStages <- as.character(preRep)
  }
  if (length(Rep) > 1) {
    repStages <- paste(Rep[1], "-", Rep[length(Rep)], sep = "")
  } else {
    repStages <- as.character(Rep)
  }
  if (length(postRep) > 1) {
    postRepStages <- paste(postRep[1], "-", postRep[length(postRep)], sep = "")
  } else {
    postRepStages <- as.character(postRep)
  }
  
  return(c(propStages, preRepStages, repStages, postRepStages))
}
