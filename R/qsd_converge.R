#' Calculate time to reach quasi-stationary stage distribution from a matrix
#' population model
#'
#' Calculate the time in the projection of a cohort through a matrix population
#' model at which a defined quasi-stationary stage distribution is reached.
#'
#' Some matrix population models are parameterised with a stasis loop at the
#' largest/most-developed stage class, which can lead to artefactual pleateaus
#' in the mortality or fertility trajectories derived from such models. These
#' plateaus occur as a projected cohort approaches its stationary stage
#' distribution (SSD). Though there is generally no single time point at which
#' the SSD is reached, we can define a quasi-stationary stage distribution (QSD)
#' based on a given distance threshold from the SSD, and calculate the number of
#' time steps required for a cohort to reach the QSD. This quantity can then be
#' used to subset age trajectories of mortality or fertility to periods earlier
#' than the QSD, so as to avoid artefactual plateaus in mortality or fertility.
#' 
#' @param mat A matrix population model, or component thereof (i.e. a square
#'   projection matrix)
#' @param start The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param conv Proportional distance threshold from the stationary stage
#'   distribution indicating convergence. E.g. This value should be 0.05 if the
#'   user wants to obtain the time step when the stage distribution is within a
#'   distance of 5\% of the stationary stage distribution.
#' @param N Maximum number of time steps over which the population will be
#'   projected. Time steps are in the same units as the matrix population model
#'   (see AnnualPeriodicity column in COM(P)ADRE). Defaults to 1000.
#' 
#' @return An integer indicating the first time step at which the
#'   quasi-stationary stage distribution is reached (or an \code{NA} and a
#'   warning if the quasi-stationary distribution is not reached).
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
#'   specifed starting stage), but often lead to non-ergodicity and therefore
#'   prevent the reliable calculation of SSD. If the reduced matrix is ergodic,
#'   the function internally updates the starting stage class and continues with
#'   the regular calculation. Otherwise, if the matrix cannot be made ergodic,
#'   the function will throw an error.
#' 
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @author Owen Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'
#'   Jones, O.R. et al. (2014) Diversity of ageing across the tree of life.
#'   Nature, 505(7482), 169–173. https://doi.org/10.1038/nature12789
#'
#'   Salguero-Gómez R. (2018) Implications of clonality for ageing research.
#'   Evolutionary Ecology, 32, 9-28. https://doi.org/10.1007/s10682-017-9923-2
#' 
#' @examples
#' mat <- rbind(c(0.1,   0,   0,   0),
#'              c(0.5, 0.2, 0.1,   0),
#'              c(  0, 0.3, 0.3, 0.1),
#'              c(  0,   0, 0.5, 0.6))
#' 
#' qsd_converge(mat)
#' 
#' @importFrom popbio stable.stage
#' @importFrom popdemo isErgodic
#' @export qsd_converge
qsd_converge <- function(mat, start = 1L, conv = 0.05, N = 1000L) {
  
  # validate arguments
  checkValidMat(mat, warn_surv_issue = TRUE)
  checkValidStartLife(start, mat)
  
  # if not ergodic, remove stages not connected from start
  if (!isErgodic(mat)) {
    
    n <- rep(0, nrow(mat))
    n[start] <- 1
    
    nonzero <- rep(FALSE, nrow(mat))
    nonzero[start] <- TRUE
    
    t <- 1L
    
    while (!all(nonzero) & t < (nrow(mat) * 3)) {
      n <- mat %*% n
      nonzero[n > 0] <- TRUE
      t <- t + 1L
    }
    mat <- as.matrix(mat[nonzero,nonzero])
    start <- which(which(nonzero) == start)
  }
  
  if (!isErgodic(mat)) {
    stop("Matrix is still non-ergodic after removing stages not connected",
         "from stage 'start'")
  }

  # stable distribution
  w <- stable.stage(mat)
  
  # set up a cohort with 1 individ in first stage class, and 0 in all others
  n <- rep(0, nrow(mat))
  n[start] <- 1
  
  # iterate cohort (n = cohort population vector, p = proportional structure)
  dist <- conv + 1
  t <- 0L

  while (!is.na(dist) & dist > conv & t < N) {
    dist <- 0.5 * (sum(abs(n - w)))
    n <- mat %*% n
    n <- n / sum(n)
    t <- t + 1L
  }
  
  return(ifelse(is.na(dist) | dist > conv, NA_integer_, t)) 
}
