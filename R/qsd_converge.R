#' Calculate time to reach quasi-stationary stage distribution
#'
#' @description
#' Calculates the time for a cohort projected with a matrix population
#' model to reach a defined quasi-stationary stage distribution.
#'
#' @param mat A matrix population model, or component thereof (i.e., a square
#'   projection matrix). Optionally with named rows and columns indicating the
#'   corresponding life stage names.
#' @param start The index (or stage name) of the first stage at which the author
#'   considers the beginning of life. Defaults to \code{1}. Alternately, a
#'   numeric vector giving the starting population vector (in which case
#'   \code{length(start)} must match \code{ncol(matU))}. See section
#'   \emph{Starting from multiple stages}.
#' @param conv Proportional distance threshold from the stationary stage
#'   distribution indicating convergence. For example, this value should be
#'   \code{0.01} if the user wants to obtain the time step when the stage
#'   distribution is within a distance of 1\% of the stationary stage
#'   distribution.
#' @param N Maximum number of time steps over which the population will be
#'   projected. Time steps are in the same units as the matrix population model
#'   (see \code{AnnualPeriodicity} column in COM(P)ADRE metadata). Defaults to
#'   \code{100,000}.
#'
#' @details
#' Some matrix population models are parameterised with a stasis loop at the
#' largest/most-developed stage class, which can lead to artefactual plateaus
#' in the mortality or fertility trajectories derived from such models. These
#' plateaus occur as a projected cohort approaches its stationary stage
#' distribution (SSD). Though there is generally no single time point at which
#' the SSD is reached, we can define a quasi-stationary stage distribution (QSD)
#' based on a given distance threshold from the SSD, and calculate the number of
#' time steps required for a cohort to reach the QSD. This quantity can then be
#' used to subset age trajectories of mortality or fertility to periods earlier
#' than the QSD, so as to avoid artefactual plateaus in mortality or fertility.
#'
#' \strong{Starting from multiple stages}
#'
#' Rather than specifying argument \code{start} as a single stage class from
#' which all individuals start life, it may sometimes be desirable to allow for
#' multiple starting stage classes. For example, if we want to start our
#' calculation of QSD from reproductive maturity (i.e., first reproduction), we
#' should account for the possibility that there may be multiple stage classes
#' in which an individual could first reproduce.
#'
#' To specify multiple starting stage classes, specify argument \code{start} as
#' the desired starting population vector, giving the proportion
#' of individuals starting in each stage class (the length of \code{start}
#' should match the number of columns in the relevant MPM).
#'
#' @seealso
#' \code{\link{mature_distrib}} for calculating the proportion of
#' individuals achieving reproductive maturity in each stage class.
#'
#'
#' @note The time required for a cohort to reach its QSD depends on the initial
#'   population vector of the cohort (for our purposes, the starting stage
#'   class), and so does not fundamentally require an ergodic matrix (where the
#'   long-term equilibrium traits are independent of the initial population
#'   vector). However, methods for efficiently calculating the stationary stage
#'   distribution (SSD) generally do require ergodicity.
#'
#'   If the supplied matrix (\code{mat}) is non-ergodic, \code{qsd_converge}
#'   first checks for stage classes with no connection (of any degree) from the
#'   starting stage class specified by argument \code{start}, and strips such
#'   stages from the matrix. These unconnected stages have no impact on
#'   age-specific traits that we might derive from the matrix (given the
#'   specified starting stage), but often lead to non-ergodicity and therefore
#'   prevent the reliable calculation of SSD. If the reduced matrix is ergodic,
#'   the function internally updates the starting stage class and continues with
#'   the regular calculation. Otherwise, if the matrix cannot be made ergodic,
#'   the function will return \code{NA} with a warning.
#'
#' @return An integer indicating the first time step at which the
#'   quasi-stationary stage distribution is reached (or an \code{NA} and a
#'   warning if the quasi-stationary distribution is not reached).
#'
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @author Owen Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#'
#' @family life tables
#'
#' @references
#' Caswell, H. 2001. Matrix Population Models: Construction, Analysis, and
#' Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968
#'
#' Horvitz, C. C., & Tuljapurkar, S. 2008. Stage dynamics, period survival, and
#' mortality plateaus. The American Naturalist, 172(2), 203–215.
#'
#' Jones, O. R., Scheuerlein, A., Salguero-Gomez, R., Camarda, C. G., Schaible,
#' R., Casper, B. B., Dahlgren, J. P., Ehrlén, J., García, M. B., Menges, E.,
#' Quintana-Ascencio, P. F., Caswell, H., Baudisch, A. & Vaupel, J. 2014.
#' Diversity of ageing across the tree of life. Nature 505, 169-173.
#' <doi:10.1038/nature12789>
#'
#' Salguero-Gomez R. 2018. Implications of clonality for ageing research.
#' Evolutionary Ecology, 32, 9-28. <doi:10.1007/s10682-017-9923-2>
#'
#' @examples
#' data(mpm1)
#'
#' # starting stage = 2 (i.e., "small")
#' qsd_converge(mpm1$matU, start = 2)
#' qsd_converge(mpm1$matU, start = "small") # equivalent using named life stages
#'
#' # convergence threshold = 0.001
#' qsd_converge(mpm1$matU, start = 2, conv = 0.001)
#'
#' # starting from first reproduction
#' repstages <- repro_stages(mpm1$matF)
#' n1 <- mature_distrib(mpm1$matU, start = 2, repro_stages = repstages)
#' qsd_converge(mpm1$matU, start = n1)
#'
#' @importFrom popdemo isErgodic project
#' @export qsd_converge
qsd_converge <- function(mat, start = 1L, conv = 0.01, N = 1e5L) {
  # validate arguments
  checkValidMat(mat)
  checkValidStartLife(start, mat, start_vec = TRUE)

  if (length(start) == 1) {
    start_vec <- rep(0.0, nrow(mat))
    if (!is.null(dimnames(mat))) {
      checkMatchingStageNames(mat)
      names(start_vec) <- colnames(mat)
    }
    start_vec[start] <- 1.0
  } else {
    start_vec <- start
  }

  # if not ergodic, remove stages not connected from start
  if (!popdemo::isErgodic(mat)) {
    nonzero <- rep(FALSE, nrow(mat))
    nonzero[start_vec > 0] <- TRUE

    n <- start_vec
    t <- 1L

    while (!all(nonzero) && t < (nrow(mat) * 2)) {
      n <- mat %*% n
      nonzero[n > 0] <- TRUE
      t <- t + 1L
    }

    mat <- as.matrix(mat[nonzero, nonzero])
    start_vec <- start_vec[nonzero]
  }

  # if still not ergodic, check whether observed dist at t = N matches stable
  #  dist
  if (!popdemo::isErgodic(mat)) {
    # check whether stable dist is __0__
    check_stable_zero <- stable_zero(mat, n1 = start_vec)
    if (check_stable_zero) {
      w <- rep(0, nrow(mat))
    } else {
      w <- popdemo::eigs(mat, what = "ss")

      n <- start_vec
      dist <- 1
      t <- 0L

      while (dist > 0.001 && t < N) {
        n <- mat %*% n
        n <- n / sum(n)
        dist <- 0.5 * (sum(abs(n - w)))
        t <- t + 1L
      }

      if (dist > 0.001) {
        warning(strwrap(prefix = " ", initial = "", "Matrix is still non-ergodic
        after removing stages not connected from stage 'start', and stable
        distribution doesn't match observed distribution after N iterations."),
          call. = FALSE
        )
        cat("\n")
        return(NA_integer_)
      }
    }
  } else {
    w <- stable.stage(mat)
  }

  # set up a cohort with 1 individ in first stage class, and 0 in all others
  n <- start_vec

  # iterate cohort (n = cohort population vector, p = proportional structure)
  dist <- conv + 1
  t <- 0L

  while (!is.na(dist) && dist > conv && t < N) {
    dist <- 0.5 * (sum(abs(n - w)))
    n <- mat %*% n
    if (sum(n) > 0) n <- n / sum(n)
    t <- t + 1L
  }

  return(ifelse(is.na(dist) | dist > conv, NA_integer_, t))
}


#' @noRd
stable_zero <- function(mat, n1) {
  n <- n1
  for (i in seq_len(nrow(mat))) {
    n <- mat %*% n
  }
  return(sum(n) == 0)
}
