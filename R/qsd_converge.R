#' Calculate time to reach quasi-stationary stage distribution from a matrix
#' population model
#'
#' Calculate the time in the projection of a cohort through a matrix population
#' model at which a defined quasi-stationary stage distribution is reached.
#'
#' Size or stage-based matrix population models (i.e. Lefkovitch models) and
#' some age-based matrix population models (i.e. Leslie models with a last class
#' greater than or equal to final age value) are typically parameterised with a
#' stasis loop in the largest/most-developed stage (e.g. adult survival). The
#' assumption of constancy in vital rates for individuals in that last stage
#' typically results in flat mortality and fertility plateaus. These plateaus
#' may result in mathematical artefacts when examining age-specific patterns
#' derived from age-from-stage matrix decompositions (Caswell 2001). The
#' Quasi-stationary Stage Distribution (QSD) can be used to circumvent this
#' problem. The QSD is the stage distribution that is reached some time before
#' the ultimate stable stage distribution (SSD, the normalised right eigenvector
#' of the transition matrix). With this approach, the user can ask what is the
#' time step in the projection at which the cohort approximates its stable stage
#' distribution with a given convergence tolerance level (e.g. 95%). This metric
#' allows the user to use only age-based information from before this point. See
#' the online supplementary information of Jones et al. (2014) for further
#' details.
#' 
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param conv Departure of the projection population vector from convergence to
#'   the stationary stage distribution. E.g. this value should be 0.05 if the
#'   user wants to obtain the time step when the stage distribution is within 5%
#'   of the stationary stable distribution.
#' @param N Number of time steps over which the population will be projected.
#'   Time steps are in the same units as the matrix population model (see
#'   AnnualPeriodicity column in COM(P)ADRE). Defaults to 1000.
#' @param ergodicFix logical: fix nonergodic survival (U) matrices, by removing 
#'   offending stages, or not?
#' @return An integer indicating the first time step at which the
#'   quasi-stationary stage distribution is reached (or an \code{NA} and a
#'   warning if the quasi-stationary distribution is not reached).
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @author Owen Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'
#'   Jones, O.R. et al. (2014) Diversity of ageing across the tree of life.
#'   Nature, 505(7482), 169–173. https://doi.org/10.1038/nature12789
#'
#'   Salguero-Gómez R. (2018) Implications of clonality for ageing research.
#'   Evolutionary Ecology, 32, 9-28. https://doi.org/10.1007/s10682-017-9923-2
#' @keywords methods
#' @examples
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#' 
#' qsd_converge(matU, 1)
#' @importFrom popbio stable.stage
#' @importFrom popdemo isErgodic
#' @export qsd_converge
qsd_converge <- function(matU, startLife = 1L, conv = 0.05, N = 1000L,
                         ergodicFix = FALSE) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(startLife, matU)
  
  # if not ergodic, remove stages not connected from startLife
  if (!isErgodic(matU) & ergodicFix == TRUE) {
    
    n <- rep(0, nrow(matU))
    n[startLife] <- 1
    
    nonzero <- rep(FALSE, nrow(matU))
    nonzero[startLife] <- TRUE
    
    t <- 1L
    
    while (!all(nonzero) & t < (nrow(matU) * 3)) {
      n <- matU %*% n
      nonzero[n > 0] <- TRUE
      t <- t + 1L
    }
    matU <- as.matrix(matU[nonzero,nonzero])
    startLife <- which(which(nonzero) == startLife)
  }
  
  if(!isErgodic(matU) & ergodicFix == FALSE){
    stop("matrix is nonergodic: cannot calculate lx from matU")
  }

  # stable distribution
  w <- stable.stage(matU)
  
  # set up a cohort with 1 individ in first stage class, and 0 in all others
  n <- rep(0, nrow(matU))
  n[startLife] <- 1
  
  # iterate cohort (n = cohort population vector, p = proportional structure)
  dist <- conv + 1
  t <- 0L

  while (!is.na(dist) & dist > conv & t < N) {
    dist <- 0.5 * (sum(abs(n - w)))
    n <- matU %*% n
    n <- n / sum(n)
    t <- t + 1L
  }
  
  return(ifelse(is.na(dist) | dist > conv, NA_integer_, t)) 
}
