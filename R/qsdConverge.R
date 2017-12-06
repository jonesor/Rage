#' Cut-off of quasi-stable distribution
#' 
#' Calculates the time in the projection of a cohort through a matrix population model at which
#' a defined quasi-stationary stage distribution is reached in \code{\link{Rage}}
#' 
#' Size or stage-based matrix population models (i.e. Lefkovitch models) and some age-based
#' matrix population models (i.e. Leslie models with a last class ≥ final age value) are
#' typically parameterised with a stasis loop in the
#' largest/most-developed stage (e.g. adult survival). The assumption
#' of constancy in vital rates for individuals in that last stage typically results in flat
#' mortality and fertility plateaus. These plateaus may result in mathematical artefacts
#' when examining age-specific patterns dervied from age-from-stage matrix
#' decompositions (Caswell 2001). The Quasi-stationary Stage Distribution (QSD)
#' can be used to circunvent this problem. The QSD is the stage distribution
#' that is reached some time before the ultimate stable stage distribution (SSD, the
#' normalised right eigenvector of the transition matrix). With this approach, the user
#' can ask what is the time step in the projectin at which the cohort
#' approximates its stable stage distribution with a given convergence tolerance level (e.g. 95%).
#' This metric allows the user to use only age-based information from before this point. See the
#' online supplementary information of Jones et al. (2014) for further details.
#' 
#' @param matU The survival-dependent matrix (a subset of the A matrix), See
#' Caswell 2001, or other references below.
#' @param conv Departure of the projection population vector from convergence
#' to the stationary st/age distribution. E.g. this value should be 0.05 if the
#' user wants to obtain the time step when the stage distribution is within 5
#' percent of the stationary st/able distribution.
#' @param startLife The stage in the life cycle of the organism, as represented
#' by the matrix population model, where the user considers 'the begining of
#' life'. This is useful to discard pre-established stages such as permanent
#' seedbanks..
#' @param nSteps Number of time steps for which the projection will be made. In
#' most cases the default of 1000 is sufficient to reach stationary st/age
#' distribution. This is on the same units as the matrix population model - see
#' MatrixPeriodicity in metadata of COMPADRE/COMADRE.
#' @return Returns a numeric value indicating the time step at which
#' quasi-stationary stage distribution has been reached.
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @author Owen Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#' Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#' 978-0878930968
#' 
#' Jones, O.R. et al. (2014) Diversity of ageing across the tree of life.
#' Nature, 505(7482), 169–173
#' 
#' Salguero-Gómez R. Implications of clonality for ageing research.
#' Evolutionary Ecology DOI 10.1007/s10682-017-9923-2
#' 
#' @keywords methods
#' @examples
#' matU <- matrix(c(0.000, 0.691, 0.000, 0.000, 0.000, 0.880, 0.000, 0.000, 0.790),nrow=3)
#' qsdConverge(matU, 0.05)
#' @export qsdConverge
#' @import popbio

qsdConverge <- function(matU, conv = 0.05, startLife = 1, nSteps = 1000){
  
  #Function to determine the cutoff age at quasi-convergence for lx and mx (Code adapted from H. Caswell's matlab code):
  
  uDim = dim(matU)
  eig = eigen.analysis(matU)
  qsd = eig$stable.stage
  qsd = as.numeric(t(matrix(qsd / sum(qsd))))
  
  #Set up a cohort
  n = rep(0, uDim[1]) #Set a population vector of zeros
  n[startLife] = 1 #Set the declared first stage of life to = 1
  
  #Iterate the cohort (n= cohort population vector, p = proportional structure)
  dist = p = NULL
  survMatrix1 <- matU
  for (j in 1:nSteps){ #j represent years of iteration
    p = n / sum(n) #Get the proportional distribution
    dist[j] = 0.5 * (sum(abs(p - qsd)))
    n = survMatrix1 %*% n #Multiply the u and n matrices to iterate
  }
  #Find the ages for convergence to conv. (default = 0.05).
  #i.e. within 5% of the QSD.
  if(min(dist, na.rm = T) < conv) {
    convage = min(which(dist < conv)) }
  if(min(dist, na.rm = T) >= conv | sum(!is.na(dist)) == 0) {
    convage = NA
    warning("Convergence not reached") }
  return(convage) 
}
