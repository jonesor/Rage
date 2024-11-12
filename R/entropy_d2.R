#' Calculate Demetrius' entropy from trajectories of age-specific survivorship
#' and fecundity
#'
#' This function calculates Demetrius' entropy from vectors of age-specific
#' survivorship (\code{lx}) and fecundity (\code{mx}). It is based on Equation
#' 4.96 of Caswell (2001).
#'
#' @section Warning: Note that this function may produce unexpected results if
#'   used on partial survivorship and fecundity trajectories. In addition, it is
#'   sensitive to the length of the these vectors. We direct users to the
#'   functions `\code{\link{shape_surv}}` and `\code{\link{shape_rep}}` which
#'   are relatively robust to these issues.
#'
#' @param lx Either a survivorship trajectory (a vector of
#'   monotonically-declining values in the interval [0,1]), or submatrix U from
#'   a matrix population model.
#' @param mx Either an age-specific fecundity trajectory (a vector of
#'   non-negative values), or submatrix U from a matrix population model.
#' @param ... Additional variables passed to `mpm_to_lx` and `mpm_to_mx` if the
#'   data are supplied as matrices. This could include the `start` argument to
#'   select a starting stage.
#'
#' @return Demetrius' entropy.
#'
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' @author Richard Hinrichsen <rich@@hinrichsenenvironmental.com>
#'
#' @family life history traits
#'
#' @references Demetrius, L., & Gundlach, V. M. 2014. Directionality theory and
#'   the entropic principle of natural selection. Entropy 16: 5428-5522.
#'   
#'   Caswell, H. 2001. Matrix Population Models: Construction, Analysis, and
#'   Interpretation. Sinauer Associates.
#'
#' @examples
#' data(mpm1)
#'
#' # derive trajectories of lx and mx, starting from stage 2
#' lx <- mpm_to_lx(mpm1$matU, start = 2)
#' mx <- mpm_to_mx(mpm1$matU, mpm1$matF, start = 2)
#'
#' # calculate Demetrius' entropy
#' entropy_d(lx, mx)
#'
#' # calculate Demetrius' entropy directly from MPM
#' entropy_d(lx = mpm1$matU, mx = mpm1$matF, start = 2)
#'
#' @export entropy_d2
entropy_d2 = function (lx, mx, ...) {
  if (inherits(lx, "matrix") && inherits(mx, "matrix")) {
    mx <- mpm_to_mx(lx, mx, ...)
  }
  if (inherits(lx, "matrix")) {
    lx <- mpm_to_lx(lx, ...)
  }
  if (any(lx < 0 | lx > 1)) {
    stop("All values of lx must be within the interval [0, 1].\n")
  }
  if (any(diff(lx) > 1e-07)) {
    stop("Values of lx must be monotonically declining.\n")
  }
  if (any(mx < 0)) {
    stop("All values of mx must be >= 0.\n")
  }
  
  lxmx <- lx * mx
  A <- get.leslie(lx, mx)
  nage <- dim(A)[1]
  lambda <- max(Mod(eigen(A)$values))
  expon <- 1:nage
  LAM <- lambda ^ (-expon)
  px2 <- lxmx * LAM
  log_px2 <- log(px2)
  log_px2[px2 == 0] <- 0
  H <- -sum(px2 * log_px2)
  
  return(H)
}