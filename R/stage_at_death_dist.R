#' Compute distribution of stage at death from a stage-structured MPM
#'
#' Calculates the expected distribution of deaths across stages for a cohort
#' governed by the survival/transition submatrix \code{U} of a matrix population
#' model. The result can be interpreted as the probability of dying in each
#' stage.
#'
#' @section Note:
#' This function uses the fundamental matrix of the Markov chain implied by
#' \code{U} to estimate the total expected time spent in each stage, then
#' combines this with per-stage death probabilities. The output describes where
#' in the life cycle deaths occur, given a specified starting stage
#' distribution.
#'
#' @section Warning:
#' Results may be misleading if \code{U} includes reproduction or if any column
#' sums exceed 1. Use only the survival/transition submatrix. The calculation
#' also assumes that all stages are transient (i.e. eventual death is certain).
#'
#' @param matU A square numeric matrix giving the survival/transition (\eqn{U})
#'   submatrix of a stage-based (Lefkovitch) matrix population model.
#' @param start Optional numeric vector of length \eqn{n} giving the initial
#'   cohort distribution across stages. It will be rescaled to sum to 1.
#'   Defaults to a uniform distribution if \code{NULL}.
#'
#' @return A numeric vector giving the expected contribution of each stage to
#'   deaths. Names are taken from \code{rownames(U)} when available.
#'
#' @author Owen Jones <jones@@biology.sdu.dk>
#'
#' @family life history traits
#'
#' @examples
#' data(mpm1)
#'
#' # Uniform starting cohort (default)
#' stage_at_death_dist(mpm1$matU)
#'
#' # Starting entirely in stage 1
#' stage_at_death_dist(mpm1$matU, start = c(1, 0, 0, 0, 0))
#'
#' # Starting with the SSD (requires the full A matrix)
#' matA <- mpm1$matU + mpm1$matF
#' ssd <- popdemo::eigs(matA, "ss")
#' stage_at_death_dist(mpm1$matU, start = ssd)
#'
#' @export
stage_at_death_dist <- function(matU, start = NULL) {
  if (!is.matrix(matU) || !is.numeric(matU) || nrow(matU) != ncol(matU)) {
    stop("matU must be a square numeric matrix")
  }

  n <- ncol(matU) # define n before using it

  if (!is.null(start)) {
    if (!is.numeric(start)) {
      stop("'start' must be numeric")
    }
    if (length(start) != n) {
      stop(paste0("'start' must have length ", n, " (number of stages)"))
    }
    if (any(start < 0)) {
      stop("'start' must not contain negative values")
    }
    if (sum(start) == 0) {
      stop("'start' must contain at least one non-zero value")
    }
  }

  w <- if (is.null(start)) rep(1 / n, n) else as.numeric(start)
  w <- w / sum(w)

  # per-stage death probability (assumes U excludes reproduction)
  d <- 1 - colSums(matU)

  # basic sanity check (optional)
  if (any(d < -1e-12)) {
    stop("Some column sums of matU exceed 1; matU should be survival/transition only.")
  }

  N <- solve(diag(n) - matU) # fundamental matrix for cohort occupancy
  occ <- as.vector(N %*% w)
  deaths <- as.vector(d * occ)

  names(deaths) <- rownames(matU)
  return(deaths)
}
