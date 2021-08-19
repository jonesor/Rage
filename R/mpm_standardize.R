#' Transform a matrix population model to a standardized form
#'
#' Transform a matrix population model to a standardized set of stage classes
#' (e.g., propagule, pre-reproductive, reproductive, and post-reproductive). The
#' transition rates in the standardized matrix are a weighted mean of the
#' transition rates and per-capita reproductive values from the relevant stages
#' of the original matrix, weighted by the relative proportion of each stage
#' class expected at the stable distribution.
#' 
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression).
#' @param matF The sexual component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions due to sexual reproduction).
#' @param matC The clonal component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions due to clonal reproduction).
#'   Defaults to \code{NULL}, indicating no clonal reproduction (i.e.
#'   \code{matC} is a matrix of zeros).
#' @param repro_stages Logical vector of length \code{ncol(matU)} indicating 
#'   which stages are reproductive. Alternatively, a vector of stage indices or 
#'   stage names of the reproductive classes.
#' @param matrix_stages Character vector of matrix stage types (e.g., "propagule",
#'   "active", or "dormant").
#' @return A list with four elements reflecting the standardized matrix and
#'   its components:
#'   \item{matA}{Standardized projection matrix}
#'   \item{matU}{Survival component of the standardized projection matrix}
#'   \item{matF}{Sexual reproduction component of the standardized projection matrix}
#'   \item{matC}{Clonal reproduction component of the standardized projection matrix}
#'   
#' @section Missing Stages:
#' The returned standardized matrix will always be of dimension \code{4}, even
#' if one or more standardized stages is missing from the original matrix
#' population model. If a standardized stage is missing, the corresponding
#' row/column of the standardized matrix will be coerced to \code{NA}.
#' 
#' @details This function is a wrapper for the functions
#'   \code{\link{mpm_rearrange}}, \code{\link{standard_stages}} and
#'   \code{\link{mpm_collapse}}, which it calls in sequence.
#' @note The method used by this function to collapse a matrix population model
#'   preserves the equilibrium population growth rate (\eqn{\lambda}) and relative
#'   stable distribution, but is not expected to preserve other demographic characteristics
#'   such as relative reproductive value, sensitivities, net reproductive rate, life
#'   expectancy, etc.
#'   
#' @author Rob Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' 
#' @family transformation
#' 
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
#' matC <- rbind(c(  0, 0.6,   0, 0.5,   0),
#'               c(  0, 0.1,   0, 0.3,   0),
#'               c(  0,   0,   0,   0,   0),
#'               c(  0,   0,   0,   0,   0),
#'               c(  0,   0,   0,   0,   0))
#' 
#' repro_stages <- c(2, 4)
#' matrix_stages <- c('prop', 'active', 'active', 'active', 'active')
#'
#' mpm_standardize(matU, matF, matC, repro_stages, matrix_stages)
#' 
#' @export mpm_standardize
mpm_standardize <- function(matU, matF, matC = NULL, repro_stages,
                            matrix_stages) {
  
  checkValidMat(matU)
  checkValidMat(matF)
  checkMatchingStageNames(matU, matF)
  if (!is.null(matC)) {
    checkValidMat(matC, warn_all_zero = FALSE)
    checkMatchingStageNames(matU, matC)
  }
  checkValidStages(matU, repro_stages)
  
  # note argument validation done by component functions
  
  # populate matC NULL, populate with zeros
  if (is.null(matC)) {
    matC <- matrix(0, nrow = ncol(matF), ncol = ncol(matF))
  }
  
  # put non-reproductive stages at the end of the matrix
  rearr <- mpm_rearrange(matU, matF, matC, repro_stages, matrix_stages)
  
  # defines which columns need to be collapsed for each of the four stages
  collapse <- standard_stages(rearr$matF,
                              rearr$repro_stages,
                              rearr$matrix_stages)
  
  # collapse
  out <- mpm_collapse(matU, matF, matC, collapse = collapse)

  return(out)
}

#' @rdname mpm_standardize
#' @export mpm_standardise
mpm_standardise <- mpm_standardize