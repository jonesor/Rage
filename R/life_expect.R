#' Calculate life expectancy from a matrix population model
#'
#' Applies Markov chain approaches to obtain mean life expectancy from a matrix
#' population model.
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param start The index of the first stage at which the author considers the
#'   beginning of life. Defaults to 1. Alternately, a numeric vector giving the
#'   starting population vector (in which case \code{length(start)} must match
#'   \code{ncol(matU))}. See section \emph{Starting from multiple stages}.
#' 
#' @return Returns life expectancy. If \code{matU} is singular (often indicating
#'   infinite life expectancy), returns \code{NA}.
#'   
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <hcaswell@@whoi.edu>
#' 
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' 
#' @section Starting from multiple stages:
#' Rather than specifying argument \code{start} as a single stage class from
#' which all individuals start life, it may sometimes be desirable to allow for
#' multiple starting stage classes. For example, if we want to start our
#' calculation of life expectancy from reproductive maturity (i.e. first
#' reproduction), we should account for the possibility that there may be
#' multiple stage classes in which an individual could first reproduce.
#' 
#' To specify multiple starting stage classes, specify argument \code{start} as
#' the desired starting population vector (\strong{n1}), giving the proportion
#' of individuals starting in each stage class (the length of \code{start}
#' should match the number of columns in the relevant MPM).
#' 
#' See function \code{\link{mature_distrib}} for calculating the proportion of
#' individuals achieving reproductive maturity in each stage class.
#' 
#' @examples
#' data(mpm1)
#' 
#' # life expectancy starting from stage class 2 
#' life_expect(mpm1$matU, start = 2)
#' 
#' # life expectancy starting from first reproduction
#' repstages <- repro_stages(mpm1$matF)
#' n1 <- mature_distrib(mpm1$matU, start = 2, repro_stages = repstages)
#' life_expect(mpm1$matU, start = n1)
#'
#' @export life_expect
life_expect <- function(matU, start = 1L) {

  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU, start_vec = TRUE)

  # matrix dimension
  matDim <- nrow(matU)
  
  if (length(start) > 1) {
    start_vec <- start / sum(start)
  } else {
    start_vec <- rep(0.0, matDim)
    start_vec[start] <- 1.0
  }
  
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)
  
  # check for errors due to singular matrix
  # if singular, return NA
  if (class(N) == "try-error" && grepl("singular", N[1])) {
    life_expect <- NA_real_
  } else {
    life_expect <- sum(colSums(N) * start_vec)
  }
  
	return(life_expect)
}

