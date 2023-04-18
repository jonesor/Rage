#' Rearrange stages of a matrix population model to segregate reproductive and
#' non-reproductive stages
#'
#' Rearrange stages of a matrix population model so that all inter-reproductive
#' stages fall in the final rows/columns of the matrix. This is a preparatory
#' step to collapsing the matrix model into a standardized set of stages (e.g.,
#' propagule, pre-reproductive, reproductive, and post-reproductive).
#'
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.,
#'   progression, stasis, and retrogression)
#' @param matF The sexual component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions due to sexual reproduction)
#' @param matC The clonal component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions due to clonal reproduction).
#'   Defaults to \code{NULL}, indicating no clonal reproduction (i.e.,
#'   \code{matC} is a matrix of zeros).
#' @param repro_stages Logical vector of length \code{ncol(matU)} indicating
#'   which stages are reproductive. Alternatively, a vector of stage indices or
#'   stage names of the reproductive classes.
#' @param matrix_stages A character vector identifying organized matrix stages.
#' @return Returns a list with 6 elements:
#' \item{matU}{Rearranged survival matrix}
#' \item{matF}{Rearranged sexual reproduction matrix}
#' \item{matC}{Rearranged clonal reproduction matrix}
#' \item{matrix_stages}{Rearranged vector of organized matrix stages}
#' \item{repro_stages}{Rearranged logical vector of reproductive stages}
#' \item{nonRepInterRep}{Numeric index for any rearranged inter-reproductive
#'  stages}
#'
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#'
#' @family transformation
#'
#' @seealso \code{\link{mpm_standardize}}
#' @examples
#' matU <- rbind(
#'   c(0.1, 0, 0, 0, 0),
#'   c(0.5, 0.2, 0.1, 0, 0),
#'   c(0, 0.3, 0.3, 0.1, 0),
#'   c(0, 0, 0.4, 0.4, 0.1),
#'   c(0, 0, 0, 0.1, 0.4)
#' )
#'
#' matF <- rbind(
#'   c(0, 1.1, 0, 1.6, 0),
#'   c(0, 0.8, 0, 0.4, 0),
#'   c(0, 0, 0, 0, 0),
#'   c(0, 0, 0, 0, 0),
#'   c(0, 0, 0, 0, 0)
#' )
#'
#' repro_stages <- c(2, 4)
#' matrix_stages <- c("prop", "active", "active", "active", "active")
#'
#' mpm_rearrange(matU, matF,
#'   repro_stages = repro_stages,
#'   matrix_stages = matrix_stages
#' )
#'
#' @export mpm_rearrange
#'
mpm_rearrange <- function(matU, matF, matC = NULL, repro_stages,
                          matrix_stages) {
  # validate arguments
  checkValidMat(matU)
  checkValidMat(matF)
  if (!is.null(matC)) checkValidMat(matC, warn_all_zero = FALSE)
  # convert repro_stages indices and names to logical vector, if needed
  if (any(is.numeric(repro_stages[!is.na(repro_stages)]))) {
    temprs <- vector(mode = "logical", length = ncol(matU))
    temprs[repro_stages] <- TRUE
    repro_stages <- temprs
  } else if (any(is.character(repro_stages[!is.na(repro_stages)]))) {
    repro_stages <- colnames(matU) %in% repro_stages
  }
  if (ncol(matU) != ncol(matF) ||
    ncol(matU) != length(repro_stages) ||
    length(repro_stages) != length(matrix_stages)) {
    stop("Arguments do not correspond to MPM of single dimension.\n",
      call. = FALSE
    )
  }
  checkMatchingStageNames(matU, matF)
  if (!is.null(matC)) {
    checkValidMat(matC, warn_all_zero = FALSE)
    checkMatchingStageNames(matU, matC)
  }
  checkValidStages(matU, repro_stages)

  # populate matC with zeros, if NULL
  if (is.null(matC)) {
    matC <- matrix(0, nrow = ncol(matF), ncol = ncol(matF))
  }

  # preliminaries
  matDim <- ncol(matF) # matrix dim
  Rep <- which(repro_stages == TRUE) # reproductive stage indices
  out <- NULL # initalize output list

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
    out$matrix_stages <- matrix_stages[rearrange]
    out$repro_stages <- repro_stages[rearrange]
  } else {
    # else, no need to rearrange
    out$matU <- matU
    out$matF <- matF
    out$matC <- matC
    out$matrix_stages <- matrix_stages
    out$repro_stages <- repro_stages
  }

  # inter-repro stages that were moved to the last column(s)
  out$nonRepInterRep <- ifelse(length(nonRepInterRep) > 0, nonRepInterRep, NA)

  return(out)
}
