#' Calculate life expectancy from a matrix population model
#'
#' Applies Markov chain approaches to obtain mean life expectancy from a matrix
#' population model.
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param startLife Index of the first stage at which the author considers the
#'   beginning of life. Defaults to 1.
#' @return Returns life expectancy. If \code{matU} is singular (often indicating
#'   infinite life expectancy), returns \code{NA}.
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <hcaswell@@whoi.edu>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' @examples
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
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
  }
  
  # matrix dimension
  matDim <- nrow(matU)

  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)
  
  # check for errors due to singular matrix
  # if singular, return NA
  if (class(N) == 'try-error' && grepl('singular', N[1])) {
    life_expect <- NA_real_
  } else {
    life_expect <- colSums(N)[startLife]
  }
  
	return(life_expect)
}
