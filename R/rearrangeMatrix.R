#' Rearrange matrix stages to segregate reproductive and non-reproductive stages
#' 
#' Rearrange matrix stages so that all inter-reproductive stages fall in the
#' final rows/columns of the matrix. This is a preparatory step to collapsing
#' the matrix model into a standardized set of stages (e.g. propagule,
#' pre-reproductive, reproductive, and post-reproductive).
#'
#' @export
#' @param matU Survival matrix
#' @param matF Fecundity matrix
#' @param matrixStages A character vector identifying the matrix stages
#' @param reproStages Logical vector identifying which stages reproductive
#' @return Returns a list with 5 elements: the rearranged survival matrix
#'   (\code{matU}), the rearranged fecundity matrix (\code{matF}), the
#'   rearranged vector of reproductive stages (\code{reproStages}), the numeric
#'   index for any rearranged inter-reproductive stages (\code{nonRepInterRep}),
#'   and the numeric index for the maximum reproductive stage in the rearranged
#'   reproductive stage vector (\code{maxRep}).
#'   
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @examples
#' matU <- rbind(c(0, 0, 0, 0, 0), c(0.1, 0.16, 0, 0, 0), c(0.2, 0.23, 0.12, 0,
#' 0), c(0, 0, 0.34, 0.53, 0), c(0, 0, 0, 0.34, 0))
#' 
#' matF <- rbind(c(0, 0, 0, 0, 0), c(0, 0.2, 0, 0.1, 0), c(0, 0.2, 0, 0.1, 0),
#' c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
#' 
#' reproStages <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
#' matrixStages <- c('prop', 'active', 'active', 'active', 'active')
#' rearrangeMatrix(matU, matF, reproStages, matrixStages)
#' 
#' 
rearrangeMatrix <- function(matU, matF, reproStages, matrixStages) {
  if (!(identical(dim(matU), dim(matF)) && identical(ncol(matF),
                                                     length(reproStages)))) {
    stop("Expecting matrices with equal dimensions", call. = FALSE)
  }
  
  Rep <- which(reproStages == TRUE)
  
  reArrange <- NULL
  matDim <- dim(matF)[1]
  
  if (length(Rep) > 0) {
    allRep <- Rep[1]:Rep[length(Rep)]
  } else {
    allRep <- integer(0)
  }
  
  ## These are stages that are inter-reproductive but are truly non-reproductive:
  nonRepInterRep <- allRep[which(!allRep %in% Rep)]
  if (length(nonRepInterRep) > 0) {
    allElseStages <- which(!1:matDim %in% nonRepInterRep)
    reArrangeStages <- c(allElseStages, nonRepInterRep)
    reArrange$matU <- matU[reArrangeStages, reArrangeStages]
    reArrange$matF <- matF[reArrangeStages, reArrangeStages]
    reArrange$reproStages <- reproStages[reArrangeStages]
    reArrange$matrixStages <- matrixStages[reArrangeStages]
  } else {
    ## No non-repro or inter-repro stages so no need to rearrange matrices
    reArrange$matU <- matU
    reArrange$matF <- matF
    reArrange$reproStages <- reproStages
    reArrange$matrixStages <- matrixStages
  }
  
  ## Stages that were moved to the end
  reArrange$nonRepInterRep <- ifelse(length(nonRepInterRep) > 0, 
                                     nonRepInterRep, NA)
  ## Max reproductive stage after rearrangement
  rearrRep <- which(reArrange$reproStages == TRUE)
  reArrange$maxRep <- ifelse(length(rearrRep) > 0, max(rearrRep), NA)
  return(reArrange)
}
