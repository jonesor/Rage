#' Calculate Keyfitz's entropy from a trajectory of age-specific survivorship
#'
#' Calculate Keyfitz's entropy from a vector of age-specific survivorship
#' (\code{lx}), or from the U submatrix of a matrix population model.
#'
#' @section Warning:
#' Note that this function may produce unexpected results if used on partial
#' survivorship trajectories. In addition, it is sensitive to the length of the
#' survivorship vector. We direct users to the function
#' `\code{\link{shape_surv}}` which is relatively robust to these issues.
#' 
#' Furthermore, de Vries et al. 2023 have shown that the way this function
#' calculates entropy is problematic for other reasons. We recommend to use
#' `\code{\link{entropy_k_age}}` or ``\code{\link{entropy_k_stage}}` as
#' alternatives, See de Vries et al. 2023 for details.
#'
#' @param lx Either a survivorship trajectory (a vector of
#'   monotonically-declining values in the interval [0,1]), or submatrix U from
#'   a matrix population model.
#' @param trapeze A logical argument indicating whether the composite trapezoid
#'   approximation should be used for approximating the definite integral.
#' @param ... Additional variables passed to `mpm_to_lx` if data are supplied as
#'   a matrix.  This could include the `start` argument to select a starting
#'   stage.
#'
#' @return Keyfitz's life table entropy.
#'
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#'
#' @references Keyfitz, N. 1977. Applied Mathematical Demography. New York:
#'   Wiley.
#'
#'   Demetrius, L., & Gundlach, V. M. 2014. Directionality theory and
#'   the entropic principle of natural selection. Entropy 16: 5428-5522.
#'
#'   de Vries, C., Bernard, C., & Salguero-Gómez, R. 2023. Discretising
#'   Keyfitz' entropy for studies of actuarial senescence and comparative
#'   demography. Methods in Ecology and Evolution, 14, 1312–1319.
#'   https://doi.org/10.1111/2041-210X.14083
#'
#' @family life history traits
#'
#' @examples
#' data(mpm1)
#'
#' # derive lx trajectory, starting from stage 2
#' lx <- mpm_to_lx(mpm1$matU, start = 2)
#'
#' # calculate Keyfitz' entropy
#' entropy_k(lx)
#'
#' # use trapezoid approximation for definite integral
#' entropy_k(lx, trapeze = TRUE)
#'
#' # calculate directly from the matrix
#' entropy_k(mpm1$matU)
#'
#' @export entropy_k
entropy_k <- function(lx, trapeze = FALSE, ...) {
  # This function is deprecated, use entropy_k2 which has two functions within
  # it suitable for age-based and stage-based matrices, selected with a type
  # argument.
  .Deprecated(msg = paste(
    "The way `entropy_k` calculates entropy is problematic.",
    "Use `entropy_k_age` or `entropy_k_stage` as alternatives.",
    "See de Vries et al. 2023 for details.",
    "After deprecation, this function will be renamed `life_elas`."
  ))

  if (inherits(lx, "numeric")) {
    # validate arguments
    if (any(lx < 0 | lx > 1)) {
      stop("All values of lx must be within the interval [0, 1].\n")
    }
    if (any(diff(lx) > 1e-7)) {
      stop("Values of lx must be monotonically declining.\n")
    }
  }

  if (inherits(lx, "matrix")) {
    lx <- mpm_to_lx(lx, ...)
  }

  # remove ages at/beyond which lx is 0 or NA
  lx <- lx[1:max(which(lx > 0))]

  # Calculate Keyfitz's entropy
  if (trapeze == TRUE) {
    x <- seq_along(lx)
    H <- -area_under_curve(x, lx * log(lx)) / area_under_curve(x, lx)
  } else {
    H <- -sum(lx * log(lx)) / sum(lx)
  }

  return(H)
}
