#' Identify stages corresponding to different parts of the reproductive life
#' cycle
#'
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
#' @param matF The sexual component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions only due to \emph{sexual}
#'   reproduction). It assumes that it has been rearranged so that
#'   non-reproductive stages are in the final rows/columns.
#' @param repro_stages Logical vector identifying which stages are reproductive.
#' @param matrix_stages (character) vector of stages, values are \code{prop}
#' (propagule), \code{active}, and \code{dorm} (dormant).
#'
#' @return A list with four elements:
#'   \item{propStages}{Position of the propagule stages}
#'   \item{preRepStages}{Position of the pre-reproductive stages}
#'   \item{repStages}{Position of the reproductive stages}
#'   \item{postRepStages}{Position of the post-reproductive stages}
#'
#' @author Rob Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#'
#' @family transformation
#'
#' @note Dormant stages are not currently handled.
#' @seealso \code{\link{mpm_standardize}}
#'
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
#' repro_stages <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
#' matrix_stages <- c("prop", "active", "active", "active", "active")
#'
#' r <- mpm_rearrange(matU, matF,
#'   repro_stages = repro_stages,
#'   matrix_stages = matrix_stages
#' )
#'
#' standard_stages(r$matF, r$repro_stages, r$matrix_stages)
#'
#' @export standard_stages
#'
standard_stages <- function(matF, repro_stages, matrix_stages) {
  # validate arguments
  checkValidMat(matF, warn_all_zero = FALSE)
  if (ncol(matF) != length(repro_stages) ||
    length(repro_stages) != length(matrix_stages)) {
    stop("Arguments do not correspond to MPM of the same dimension.\n",
      call. = FALSE
    )
  }
  if (!any(repro_stages == TRUE)) {
    stop(strwrap(prefix = " ", initial = "", "Cannot identify standardised
    stages because no stages are reproductive (i.e., at least one element of
                 repro_stages must be TRUE).\n"),
      call. = FALSE
    )
  }

  if ("prop" %in% matrix_stages) {
    propStage <- which(matrix_stages == "prop")
  } else {
    propStage <- NA
  }

  # set max reproductive stage
  maxRep <- max(which(repro_stages == TRUE))

  # prerep
  matDim <- nrow(matF)
  Rep <- which(repro_stages == TRUE)

  if (min(Rep) == 1) {
    preRep <- NA
  } else if (!is.na(propStage[1]) && (min(Rep) - max(propStage) == 1)) {
    preRep <- NA
  } else {
    preRep <- min(which(matrix_stages == "active")):(min(Rep) - 1)
  }

  # postrep
  if (maxRep == matDim) {
    postRep <- NA
  } else {
    postRep <- (maxRep + 1):matDim
  }

  return(list(
    propStages = propStage,
    preRepStages = preRep,
    repStages = Rep,
    postRepStages = postRep
  ))
}
