#' Calculate measures of life expectancy and longevity
#'
#' A function to calculate life expectancy and maximum longevity of individuals
#' in a matrix population model
#'
#' This function applies Markov chain approaches to obtain the mean and variance
#' in life expectancy, and/or a loop to calculate maximum longevity.
#'
#' @param matU A matrix containing only survival-dependent processes (growth,
#'   stasis, shrinkage).
#' @param startLife Index of the first stage at which the author considers the
#'   beginning of life.
#' @param initPop Initial population size.
#' @param run Number of iterations that maximum longevity will be calculated
#'   over.
#' @param method Should function return life expectancy (\code{"life"}), maximum
#'   longevity (\code{"longevity"}), or both (\code{"both"})?
#' @return Returns a list including the following named elements (each of which
#'   contain a single numeric value):
#'   
#'   If \code{method == "life"}:
#'
#'   - \code{eta}: mean life expectancy conditional on entering the life cycle
#'   at life stage described by \code{startLife}
#'
#'   - \code{var_eta}: variance in life expectancy conditional on entering the
#'   life cycle at life stage described by \code{startLife}
#'   
#'   If \code{method == "longevity"}:
#'
#'   - \code{Max}: maximum longevity observed by iterating a population vector
#'   with \code{initPop} individuals in stage \code{startLife} up to \code{run}
#'   times
#'   
#'   If \code{method == "both"}, returns a list containing all three elements
#'   described above.
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <hcaswell@@whoi.edu>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'
#'   Morris, W. F., and D. F. Doak. 2003. Quantitative Conservation Biology:
#'   Theory and Practice of Population Viability Analysis. Sinauer Associates,
#'   Sunderland, Massachusetts, USA
#' @keywords ~kwd1 ~kwd2
#' @examples
#' matU <- matrix(c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1),
#'                nrow = 4, byrow = TRUE)
#' 
#' longevity(matU, startLife = 1, initPop = 100, run = 1000)
#' 
#' @export longevity
longevity <- function(matU, startLife = 1, initPop = 100, run = 1000,
                      method = 'both') {
  #Function to calculate mean life expectancy and maximum longevity from
  # H. Caswell's matlab code, and Morris & Doak:

  if (missing(matU)) {
    stop('matU missing')
  } else if (any(is.na(matU))) {
    stop('matU contains missing values')
  } else if (!method %in% c('life', 'longevity', 'both')) {
    stop("method must be one of 'life', 'longevity', or 'both'")
  } else if (all(matU == 0)) {
    warning('all elements of matU are zero')
    return(list(eta = 0, var_eta = 0, Max = 0))
  }
  
  out <- NULL
  
  matDim <- dim(matU)[1]
  
  if(method %in% c('life', 'both')) {
    #Mean and variance of life expectancy
    N <- solve(diag(matDim) - matU)
    out$eta <- colSums(N)[startLife]
    out$var_eta <- (colSums(2 * N %*% N - N) - (colSums(N) * colSums(N)))[startLife]
  }
  
  if(method %in% c('longevity', 'both')) {
    #Maximum longevity up to 'run' iterations
    popVector <- c(rep(0, matDim))
    popVector[startLife] <- initPop
    lifespanLeftover <- matrix(0, run, 1)
    for (n in 1:run)	{
      lifespanLeftover[n] <- sum(popVector)
      popVector <- matU %*% popVector
    }
    
    out$Max <- max(which(lifespanLeftover > 1))
    if (out$Max == Inf) {out$Max <- run} # don't think this condition can ever be met
  }
  
	return(out)
}
