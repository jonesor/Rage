#' Perturbation analysis of transition types within a matrix population model
#'
#' @description Calculates the summed sensitivities or elasticities for various
#' transition types within a matrix population model (MPM), including stasis,
#' retrogression, progression, fecundity, and clonality.
#'
#' Sensitivities or elasticities are calculated by perturbing elements of the
#' MPM and measuring the response of the per-capita population growth rate at
#' equilibrium (\eqn{\lambda}), or, with a user-supplied function, any other
#' demographic statistic.
#'
#' @param matU The survival component submatrix of a MPM (i.e., a square
#'   projection matrix reflecting survival-related transitions; e.g.,
#'   progression, stasis, and retrogression).
#' @param matF The sexual component submatrix of a MPM (i.e., a square
#'   projection matrix reflecting transitions due to sexual reproduction).
#' @param matC The clonal component submatrix of a MPM (i.e., a square
#'   projection matrix reflecting transitions due to clonal reproduction).
#'   Defaults to \code{NULL}, indicating no clonal reproduction possible.
#' @param posU A logical matrix of the same dimension as \code{matU}, with
#'   elements indicating whether a given \code{matU} transition is possible
#'   (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matU > 0} (see
#'   Details).
#' @param posF A logical matrix of the same dimension as \code{matF}, with
#'   elements indicating whether a given \code{matF} transition is possible
#'   (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matF > 0} (see
#'   Details).
#' @param posC A logical matrix of the same dimension as \code{matC}, with
#'   elements indicating whether a given \code{matC} transition is possible
#'   (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matC > 0} (see
#'   Details).
#' @param exclude_row A vector of row indices or stages names indicating stages
#'   for which transitions \emph{to} the stage should be excluded from
#'   perturbation analysis. Alternatively, a logical vector of length
#'   \code{nrow(matU)} indicating which stages to include \code{TRUE} or exclude
#'   \code{FALSE} from the calculation. See section \emph{Excluding stages}.
#' @param exclude_col A vector of column indices or stages names indicating
#'   stages for which transitions \emph{to} the stage should be excluded from
#'   perturbation analysis. Alternatively, a logical vector of length
#'   \code{ncol(matU)} indicating which stages to include \code{TRUE} or exclude
#'   \code{FALSE} from the calculation. See section \emph{Excluding stages}.
#' @param pert The magnitude of the perturbation (defaults to \code{1e-6}).
#' @param type An argument defining whether to return `sensitivity` or
#'   `elasticity` values. Defaults to `sensitivity`.
#' @param demog_stat An argument defining which demographic statistic should be
#'   used, as in "the sensitivity/elasticity of \code{demog_stat} to matrix
#'   element perturbations." Defaults to the per-capita population growth rate
#'   at equilibrium (\eqn{lambda}). Also accepts a user-supplied function that
#'   performs a calculation on a MPM and returns a single numeric value.
#' @param ... Additional arguments passed to the function \code{demog_stat}.
#'
#' @details A transition rate of \code{0} within a matrix population model can
#' either indicate that the transition is not possible in the given life cycle
#' (e.g., tadpoles never revert to eggs), or that the transition is possible but
#' was estimated to be \code{0} in the relevant population and time period.
#' Because transition rates of zero \emph{do} generally yield non-zero
#' sensitivities, it is important to distinguish between structural (i.e.
#' impossible) zeros and sampled zeros when summing multiple sensitivities for a
#' given process (e.g., progression/growth).
#'
#' By default, the \code{perturb_} functions assume that a transition rate of
#' \code{0} indicates an impossible transition, in which case the sensitivity
#' for that transition will not be included in any calculation. Specifically,
#' the arguments \code{posX} are specified by the logical expression \code{(matX
#' > 0)}. If the matrix population model includes transitions that are possible
#' but estimated to be \code{0}, users should specify the \code{posX}
#' argument(s) manually.
#'
#' If there are no possible transitions for a given process (e.g., clonality, in
#' many species), the value of sensitivity or elasticity returned for that
#' process will be \code{NA}.
#'
#' @return A list with 5 elements: \item{stasis}{The sensitivity or elasticity
#'   of \code{demog_stat} to stasis.} \item{retrogression}{The sensitivity or
#'   elasticity of \code{demog_stat} to retrogression.} \item{progression}{The
#'   sensitivity or elasticity of \code{demog_stat} to progression.}
#'   \item{fecundity}{The sensitivity or elasticity of \code{demog_stat} to
#'   sexual fecundity.} \item{clonality}{The sensitivity or elasticity of
#'   \code{demog_stat} to clonality.}
#'
#' @section Excluding stages: It may be desirable to exclude one or more stages
#'   from the calculation. For instance, we might not believe that 'progression'
#'   to a dormant stage class truly reflects progression. In this case we could
#'   exclude transitions \emph{to} the dormant stage class using the argument
#'   \code{exclude_row}. We may or may not want to ignore progression
#'   transitions \emph{from} the dormant stage class, which can be done in a
#'   similar way using the argument \code{exclude_col}. The \code{exclude_}
#'   arguments simply set the relevant row or column of the \code{posX}
#'   arguments to \code{FALSE}, to prevent those transitions from being used in
#'   subsequent calculations.
#'
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#'
#' @family perturbation analysis
#'
#' @examples
#' matU <- rbind(
#'   c(0.1, 0, 0, 0),
#'   c(0.5, 0.2, 0.1, 0),
#'   c(0, 0.3, 0.3, 0.1),
#'   c(0, 0, 0.5, 0.6)
#' )
#'
#' matF <- rbind(
#'   c(0, 0, 1.1, 1.6),
#'   c(0, 0, 0.8, 0.4),
#'   c(0, 0, 0, 0),
#'   c(0, 0, 0, 0)
#' )
#'
#'
#' perturb_trans(matU, matF)
#'
#' # Use a larger perturbation than the default of 1e-6.
#' perturb_trans(matU, matF, pert = 0.01)
#'
#' # Calculate the sensitivity/elasticity of the damping ratio to perturbations.
#' # First, define function for damping ratio:
#' damping <- function(matA) {
#'   eig <- eigen(matA)$values
#'   dm <- rle(Mod(eig))$values
#'   return(dm[1] / dm[2])
#' }
#'
#' # Second, run the perturbation analysis using demog_stat = "damping".
#' perturb_trans(matU, matF, demog_stat = "damping")
#'
#' @export perturb_trans
perturb_trans <- function(matU, matF, matC = NULL,
                          posU = matU > 0, posF = matF > 0, posC = matC > 0,
                          exclude_row = NULL, exclude_col = NULL,
                          pert = 1e-6, type = "sensitivity",
                          demog_stat = "lambda", ...) {
  # Validate arguments
  checkValidMat(matU)
  checkValidMat(matF)

  if (!is.null(matC)) {
    checkValidMat(matC, warn_all_zero = FALSE)
  }
  checkValidStages(matU, exclude_row)
  checkValidStages(matU, exclude_col)
  type <- match.arg(type, c("sensitivity", "elasticity"))

  # Get statfun
  if (is.character(demog_stat) && demog_stat == "lambda") {
    statfun <- lambda
  } else {
    statfun <- try(match.fun(demog_stat), silent = TRUE)
    if (inherits(statfun, "try-error")) {
      stop(strwrap(
        prefix = " ", initial = "",
        "`demog_stat` must be `lambda` or the name
      of a function that returns a single numeric value.\n"
      ), call. = FALSE)
    }
  }

  # Matrix dimension
  m <- nrow(matU)

  # If matC null, convert to zeros
  if (is.null(matC)) {
    matC <- matrix(0, m, m)
    posC <- matrix(FALSE, m, m)

    if (!is.null(dimnames(matU))) {
      dimnames(posC) <- dimnames(matC) <- dimnames(matU)
    }
  }

  # Excluded stage classes
  posU[exclude_row, ] <- FALSE
  posU[, exclude_col] <- FALSE
  posF[exclude_row, ] <- FALSE
  posF[, exclude_col] <- FALSE
  posC[exclude_row, ] <- FALSE
  posC[, exclude_col] <- FALSE

  # Combine components into matA
  matA <- matU + matF + matC

  # Lower and upper triangles (reflecting growth and retrogression)
  lwr <- upr <- matrix(FALSE, nrow = m, ncol = m)
  lwr[lower.tri(lwr)] <- TRUE
  upr[upper.tri(upr)] <- TRUE

  posStasi <- posU & diag(m)
  posRetro <- posU & upr
  posProgr <- posU & lwr

  # Get sensitivity or elasticity
  pertMat <- perturb_matrix(
    matA = matA, pert = pert, type = type,
    demog_stat = statfun, ...
  )

  if (type == "sensitivity") {
    stasis <- ifelse(!any(posStasi), NA_real_, sum(pertMat[posStasi]))
    retro <- ifelse(!any(posRetro), NA_real_, sum(pertMat[posRetro]))
    progr <- ifelse(!any(posProgr), NA_real_, sum(pertMat[posProgr]))
    fecund <- ifelse(!any(posF), NA_real_, sum(pertMat[posF]))
    clonal <- ifelse(!any(posC), NA_real_, sum(pertMat[posC]))
  } else {
    propU <- matU / matA
    propU[!posU] <- NA_real_
    propU[matA == 0 & posU] <- 1

    propProgr <- propRetro <- propU
    propProgr[upper.tri(propU, diag = TRUE)] <- NA
    propRetro[lower.tri(propU, diag = TRUE)] <- NA

    propStasi <- matrix(NA_real_, nrow = m, ncol = m)
    diag(propStasi) <- diag(propU)

    propF <- matF / matA
    propF[!posF] <- NA_real_
    propF[matA == 0 & posF] <- 1

    propC <- matC / matA
    propC[!posC] <- NA_real_
    propC[matA == 0 & posC] <- 1

    stasis <- sum_elast(pertMat, posStasi, propStasi)
    retro <- sum_elast(pertMat, posRetro, propRetro)
    progr <- sum_elast(pertMat, posProgr, propProgr)
    fecund <- sum_elast(pertMat, posF, propF)
    clonal <- sum_elast(pertMat, posC, propC)
  }

  return(list(
    stasis = stasis, retro = retro, progr = progr,
    fecundity = fecund, clonality = clonal
  ))
}



#  A convenience function to sum elasticities given the perturbation matrix, the
#  matrix of possible transitions, and the matrix reflecting the proportional
#  contribution of the given process to the given element.
sum_elast <- function(pert_mat, pos_mat, prop_mat) {
  ifelse(!any(pos_mat), NA_real_, sum(pert_mat * prop_mat, na.rm = TRUE))
}
