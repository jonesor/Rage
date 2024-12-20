#' Calculate Demetrius' entropy from trajectories of age-specific survivorship
#' and fecundity
#'
#' This function calculates Demetrius' entropy from vectors of age-specific
#' survivorship (\code{lx}) and fecundity (\code{mx}). Users can choose between
#' the scaled (Caswell, 2001 eqns. 4.94-4.97) or unscaled (from Demetrius 1978)
#' method.
#' 
#' The scaled version accounts for population growth or shrinkage by adjusting
#' the contributions of survivorship and fecundity using the dominant eigenvalue
#' (lambda). Specifically, each contribution is weighted by lambda raised
#' to the negative power of age. Conversely, the unscaled version does not
#' account for population growth. It calculates entropy directly from the
#' proportional contributions of survivorship and fecundity without adjustment
#' for population dynamics.
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
#'   non-negative values), or submatrix F from a matrix population model.
#' @param type Calculation type, either `scaled` (default) or `unscaled`.
#' @param ... Additional variables passed to `mpm_to_lx` and `mpm_to_mx` if the
#'   data are supplied as matrices. This could include the `start` argument to
#'   select a starting stage.
#'
#' @return Demetrius' entropy.
#'
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' @author Richard Hinrichsen <rich@@hinrichsenenvironmental.com>
#' @author Owen Jones <jones@@biology.sdu.dk>
#' 
#' @family life history traits
#'
#' @references Demetrius, L. 1978. Adaptive value, entropy and survivorship
#'   curves. Nature, 275(5677), 213–214.
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
#' entropy_d(lx, mx, type = "unscaled")
#' entropy_d(lx, mx, type = "scaled")
#'
#'
#' # calculate entropy directly from MPM
#' entropy_d(lx = mpm1$matU, mx = mpm1$matF, start = 2)
#' 
#' @export
entropy_d <- function(lx, mx, type = "scaled", ...) {
  if (inherits(lx, "matrix") && inherits(mx, "matrix")) {
    mx <- mpm_to_mx(lx, mx, ...)
  }
  
  if (inherits(lx, "matrix")) {
    lx <- mpm_to_lx(lx, ...)
  }
  
  # Validate arguments
  if (any(lx < 0 | lx > 1)) {
    stop("All values of lx must be within the interval [0, 1].\n")
  }
  if (any(diff(lx) > 1e-7)) {
    stop("Values of lx must be monotonically declining.\n")
  }
  if (any(mx < 0)) {
    stop("All values of mx must be >= 0.\n")
  }
  
  lxmx <- lx * mx
  
  if (type == "scaled") {
    # Create Leslie matrix and calculate lambda
    A <- get.leslie(lx, mx)
    nage <- dim(A)[1]
    lambda <- max(Mod(eigen(A)$values))
    expon <- 1:nage
    LAM <- lambda ^ (-expon)
    px <- lxmx * LAM
  } else if (type == "unscaled") {
    px <- lxmx / sum(lxmx)
  } else {
    stop("Invalid type. Please choose either 'scaled' or 'unscaled'.")
  }
  
  log_px <- log(px)
  log_px[px == 0] <- 0
  H <- -sum(px * log_px)
  
  return(H)
}
