#' Calculate generation time from a matrix population model
#'
#' @description
#' Calculate generation time from a matrix population model. Multiple
#' definitions of the generation time are supported: the time required for a
#' population to increase by a factor of R0 (the net reproductive rate; Caswell
#' (2001), section 5.3.5), the average parent-offspring age difference (Bienvenu
#' & Legendre (2015)), or the expected age at reproduction for a cohort (Coale
#' (1972), p. 18-19).
#'
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.,
#'   progression, stasis, and retrogression).
#' @param matR The reproductive component of a matrix population model (i.e., a
#'   square projection matrix only reflecting transitions due to reproduction;
#'   either sexual, clonal, or both).
#' @param method The method used to calculate generation time. Defaults to "R0".
#'   See Details for explanation of calculations.
#' @param ... Additional arguments passed to \code{net_repro_rate} when
#'   \code{method = "R0"} or \code{mpm_to_*} when \code{method = "cohort"}.
#'   Ignored when \code{method = "age_diff"}
#'
#' @details
#' There are multiple definitions of generation time, three of which are
#' implemented by this function:
#'
#' 1. \code{"R0"} (default): This is the number of time steps required for the
#' population to grow by a factor of its net reproductive rate, equal to
#' \code{log(R0) / log(lambda)}. Here, \code{R0} is the net reproductive rate
#' (the per-generation population growth rate; Caswell 2001, Sec. 5.3.4), and
#' \code{lambda} is the population growth rate per unit time (the dominant
#' eigenvalue of \code{matU + matR}).
#'
#' 2. \code{"age_diff"}: This is the average age difference between parents and
#' offspring, equal to \code{(lambda v w) / (v matR w)} (Bienvenu & Legendre
#' (2015)). Here, \code{lambda} is the population growth rate per unit time (the
#' dominant eigenvalue of \code{matU + matR}), \code{v} is a row vector of
#' stage-specific reproductive values (the left eigenvector corresponding to
#' \code{lambda}), and \code{w} is a column vector of the stable stage
#' distribution (the right eigenvector corresponding to \code{lambda}).
#'
#' 3. \code{"cohort"}: This is the age at which members of a cohort are expected
#' to reproduce, equal to \code{sum(x lx mx) / sum(lx mx)} (Coale (1972), p.
#' 18-19). Here, \code{x} is age, \code{lx} is age-specific survivorship, and
#' \code{mx} is age-specific fertility. See functions \code{mpm_to_lx} and
#' \code{mpm_to_mx} for details about the conversion of matrix population models
#' to life tables.
#'
#' @return Returns generation time. If \code{matU} is singular (often indicating
#'   infinite life expectancy), returns \code{NA}.
#'
#' @note Note that the units of time in returned values are the same as the
#'   projection interval (`ProjectionInterval`) of the MPM.
#'
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' @author William Petry <wpetry@@ncsu.edu>
#'
#' @family life history traits
#'
#' @references Bienvenu, F. & Legendre, S. 2015. A New Approach to the
#'   Generation Time in Matrix Population Models. The American Naturalist 185
#'   (6): 834–843. doi:10.1086/681104.
#'
#' Caswell, H. 2001. Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'
#' Coale, A.J. 1972. The Growth and Structure of Human Populations. Princeton
#' University Press. ISBN: 978-0691093574
#'
#' @examples
#' data(mpm1)
#'
#' # calculate generation time
#' gen_time(matU = mpm1$matU, matR = mpm1$matF) # defaults to "R0" method
#' gen_time(matU = mpm1$matU, matR = mpm1$matF, method = "age_diff")
#' gen_time(
#'   matU = mpm1$matU, matR = mpm1$matF, method = "cohort", lx_crit =
#'     0.001
#' )
#'
#' @export gen_time
gen_time <- function(matU, matR, method = c("R0", "age_diff", "cohort"), ...) {
  method <- match.arg(method)
  # leave remaining arg validation to functions that depend on them

  if (method == "R0") {
    R0 <- net_repro_rate(matU, matR, ...)
    lam <- lambda(matU + matR)
    out <- log(R0) / log(lam)
  } else if (method == "age_diff") {
    ignore <- list(...)
    # check for singularity
    matDim <- nrow(matU)
    N <- try(solve(diag(matDim) - matU), silent = TRUE)
    if (inherits(N, "try-error") && grepl("singular", N[1], fixed = TRUE)) {
      out <- NA_real_
    } else {
      A <- matU + matR
      lam <- lambda(A)
      v <- reproductive.value(A)
      w <- stable.stage(A)
      out <- as.numeric((lam * v %*% w) / (v %*% matR %*% w))
    }
  } else if (method == "cohort") {
    # check for singularity
    matDim <- nrow(matU)
    N <- try(solve(diag(matDim) - matU), silent = TRUE)
    if (inherits(N, "try-error") && grepl("singular", N[1], fixed = TRUE)) {
      out <- NA_real_
    } else {
      lx <- mpm_to_lx(matU, ...)
      mx <- mpm_to_mx(matU, matR, ...)
      x <- seq_along(lx)
      out <- sum(x * lx * mx) / sum(lx * mx)
    }
  } else {
    stop("Unsupported method. Must be either 'R0', 'age_diff', or 'cohort'.")
  }
  return(out)
}
