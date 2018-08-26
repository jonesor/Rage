#' Calculate life expectancy from a matrix population model
#'
#' Applies Markov chain approaches to obtain mean life expectancy from a matrix
#' population model
#'
#' @param matU A matrix containing only survival-dependent processes (growth,
#'   stasis, shrinkage).
#' @param startLife Index of the first stage at which the author considers the
#'   beginning of life. Defaults to 1.
#' @return Returns life expectancy.
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <hcaswell@@whoi.edu>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'
#'   Morris, W. F., and D. F. Doak. 2003. Quantitative Conservation Biology:
#'   Theory and Practice of Population Viability Analysis. Sinauer Associates,
#'   Sunderland, Massachusetts, USA
#' @examples
#' matU <- rbind(c(0.0, 0.0, 0.0, 0.0),
#'               c(0.5, 0.0, 0.0, 0.0),
#'               c(0.0, 0.3, 0.0, 0.0),
#'               c(0.0, 0.0, 0.1, 0.1))
#'
#' lifeExpectancy(matU)
#' lifeExpectancy(matU, startLife = 2)
#'
#' @export lifeExpectancy
lifeExpectancy <- function(matU, startLife = 1) {

  # input validation
  if (any(is.na(matU))) {
    stop('matU contains missing values')
  }
  if (all(matU == 0)) {
    warning('all elements of matU are zero')
    return(0)
  }
  
  # calculate life expectancy
  matDim <- nrow(matU)
  N <- solve(diag(matDim) - matU)
  life_expect <- colSums(N)[startLife]
  
	return(life_expect)
}
