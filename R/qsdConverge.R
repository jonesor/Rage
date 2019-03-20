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
#' @param conv Departure of the projection population vector from convergence to
#'   the stationary stage distribution. E.g. this value should be 0.05 if the
#'   user wants to obtain the time step when the stage distribution is within 5%
#'   of the stationary stable distribution.
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param N Number of time steps over which the population will be projected.
#'   Time steps are in the same units as the matrix population model (see
#'   AnnualPeriodicity column in COM(P)ADRE). Defaults to 1000.
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
#' qsdConverge(matU, 0.05)
#' @importFrom popbio stable.stage
#' @export qsdConverge
qsdConverge <- function(matU, conv = 0.05, startLife = 1, N = 1000) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(startLife, matU)
  
  # stable distribution
  w <- stable.stage(matU)
  
  # set up a cohort with 1 individ in first stage class, and 0 in all others
  n <- rep(0, nrow(matU))
  n[startLife] <- 1
  
  # iterate cohort (n = cohort population vector, p = proportional structure)
  dist <- NULL
  
  for (j in 1:N) {            # j represent years of iteration
    p <- n / sum(n)           # get the proportional distribution
    dist[j] <- 0.5 * (sum(abs(p - w)))  # distance to the stable distribution
    n <- matU %*% n           # multiply matU %*% n to iterate
  }
  
  # find time to convergence (default conv = 0.05; i.e. within 5% of w)
  if (min(dist, na.rm = TRUE) < conv) {
    convage <- min(which(dist < conv))
  }
  if (min(dist, na.rm = TRUE) >= conv) {
    convage <- NA_integer_
    warning("Convergence not reached within N")
  }
  
  return(convage) 
}
