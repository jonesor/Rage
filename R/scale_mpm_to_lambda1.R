#' Scale a matrix population model so that lambda equals 1
#'
#' Scales a matrix population model (MPM) so that the dominant eigenvalue
#' (\eqn{\lambda}) of the full projection matrix equals 1. The function accepts
#' either a full projection matrix (\code{matA}) or a set of component matrices
#' consisting of the survival component (\code{matU}) plus at least one
#' reproductive component (\code{matF} and/or \code{matC}).
#'
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.,
#'   progression, stasis, and retrogression). Required when using component
#'   inputs.
#' @param matF The sexual reproduction component of a matrix population model.
#'   Optional if \code{matC} is supplied. If omitted while \code{matC} is
#'   supplied, it is assumed to be a zero matrix.
#' @param matC The clonal reproduction component of a matrix population model.
#'   Optional if \code{matF} is supplied. If omitted while \code{matF} is
#'   supplied, it is assumed to be a zero matrix.
#' @param matA The full projection matrix of a matrix population model.
#'   Optional if component matrices are supplied. If supplied together with
#'   components, it must equal \code{matU + matF + matC}.
#'
#' @return A list with six elements:
#' \describe{
#'   \item{matA}{The scaled full projection matrix.}
#'   \item{matU}{The scaled survival matrix. Returned as \code{NULL} if only
#'   \code{matA} is supplied.}
#'   \item{matF}{The scaled sexual reproduction matrix. Returned as \code{NULL}
#'   if only \code{matA} is supplied.}
#'   \item{matC}{The scaled clonal reproduction matrix. Returned as \code{NULL}
#'   if only \code{matA} is supplied.}
#'   \item{lambda_original}{The dominant eigenvalue of the original projection
#'   matrix.}
#'   \item{lambda_scaled}{The dominant eigenvalue of the scaled projection
#'   matrix.}
#' }
#'
#' @details When component matrices are used, \code{matU} must be supplied
#' together with at least one of \code{matF} or \code{matC}. Missing
#' reproductive components are treated as zero matrices of the appropriate
#' dimension.
#'
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#'
#' @family transformation
#'
#' @examples
#' data(mpm1)
#'
#' # scale from the full projection matrix
#' matA <- mpm1$matU + mpm1$matF
#' scale_mpm_to_lambda1(matA = matA)
#'
#' # scale from component matrices, assuming no clonal reproduction
#' scale_mpm_to_lambda1(matU = mpm1$matU, matF = mpm1$matF)
#'
#' @export
scale_mpm_to_lambda1 <- function(matU = NULL,
                                 matF = NULL,
                                 matC = NULL,
                                 matA = NULL) {
  has_A <- !is.null(matA)
  has_U <- !is.null(matU)
  has_F <- !is.null(matF)
  has_C <- !is.null(matC)
  has_components <- has_U || has_F || has_C

  if (!has_A && !(has_U && (has_F || has_C))) {
    stop(
      "Supply either 'matA' or 'matU' plus at least one of 'matF' or 'matC'.",
      call. = FALSE
    )
  }

  if (has_components && !has_U) {
    stop("Argument 'matU' is required when supplying component matrices.",
      call. = FALSE
    )
  }

  validate_scale_input <- function(M, arg) {
    if (!is.matrix(M) || !is.numeric(M) || nrow(M) != ncol(M)) {
      stop("Argument '", arg, "' must be a square numeric matrix.",
        call. = FALSE
      )
    }
    if (anyNA(M)) {
      stop("Argument '", arg, "' contains missing values (i.e. <NA>).",
        call. = FALSE
      )
    }
    if (!all(is.finite(M))) {
      stop("Argument '", arg, "' must contain only finite values.",
        call. = FALSE
      )
    }
  }

  if (has_components) {
    validate_scale_input(matU, "matU")

    if (!has_F && !has_C) {
      stop(
        "Supply at least one of 'matF' or 'matC' when using component inputs.",
        call. = FALSE
      )
    }

    if (has_F) {
      validate_scale_input(matF, "matF")
      if (!all(dim(matU) == dim(matF))) {
        stop(
          "Arguments 'matU' and 'matF' must have the same dimensions.",
          call. = FALSE
        )
      }
      checkMatchingStageNames(matU, matF)
    } else {
      matF <- matrix(0, nrow = nrow(matU), ncol = ncol(matU))
      dimnames(matF) <- dimnames(matU)
    }

    if (has_C) {
      validate_scale_input(matC, "matC")
      if (!all(dim(matU) == dim(matC))) {
        stop(
          "Arguments 'matU' and 'matC' must have the same dimensions.",
          call. = FALSE
        )
      }
      checkMatchingStageNames(matU, matC)
    } else {
      matC <- matrix(0, nrow = nrow(matU), ncol = ncol(matU))
      dimnames(matC) <- dimnames(matU)
    }

    matA_from_components <- matU + matF + matC
  }

  if (has_A) {
    validate_scale_input(matA, "matA")
    if (has_components) {
      if (!all(dim(matA) == dim(matU))) {
        stop(
          "Arguments 'matA' and 'matU' must have the same dimensions.",
          call. = FALSE
        )
      }
      checkMatchingStageNames(matU, matA)
      if (!isTRUE(all.equal(matA, matA_from_components,
        tolerance = sqrt(.Machine$double.eps)
      ))) {
        stop(
          "Argument 'matA' must equal 'matU + matF + matC'.",
          call. = FALSE
        )
      }
    }
  } else {
    matA <- matA_from_components
  }

  lambda_original <- suppressWarnings(as.numeric(Re(lambda(matA))))

  if (!is.finite(lambda_original) || lambda_original <= 0) {
    stop(
      "Argument 'matA' must have a finite positive dominant eigenvalue.",
      call. = FALSE
    )
  }

  matA_scaled <- matA / lambda_original

  out <- list(
    matA = matA_scaled,
    matU = NULL,
    matF = NULL,
    matC = NULL,
    lambda_original = lambda_original,
    lambda_scaled = suppressWarnings(as.numeric(Re(lambda(matA_scaled))))
  )

  if (has_components) {
    out$matU <- matU / lambda_original
    out$matF <- matF / lambda_original
    out$matC <- matC / lambda_original
  }

  return(out)
}
