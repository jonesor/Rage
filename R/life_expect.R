#' Calculate mean and variance of life expectancy from a matrix population model
#'
#' Applies Markov chain approaches to obtain mean and variance of life
#' expectancy from a matrix population model (MPM).
#'
#' @param matU The survival component of a MPM (i.e., a square projection matrix
#'   reflecting survival-related transitions; e.g., progression, stasis, and
#'   retrogression). Optionally with named rows and columns indicating the
#'   corresponding life stage names.
#' @param start The index (or stage name) of the first stage of the life cycle
#'   which the user considers to be the beginning of life. Defaults to \code{1}.
#'   Alternately, a numeric vector giving the starting population vector (in which
#'    case \code{length(start)} must match \code{ncol(matU))}. See section
#'   \emph{Starting from multiple stages}.
#' 
#' @return Returns life expectancy. If \code{matU} is singular (often indicating
#'   infinite life expectancy), returns \code{NA}.
#'   
#' @note Note that the units of time in returned values are the same as the
#'   projection interval (`ProjectionInterval`) of the MPM.
#'   
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <hcaswell@@whoi.edu>
#' 
#' @family life history traits
#' 
#' @references Caswell, H. 2001. Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' 
#' @section Starting from multiple stages:
#' Rather than specifying argument \code{start} as a single stage class from
#' which all individuals start life, it may sometimes be desirable to allow for
#' multiple starting stage classes. For example, if the user wants to start their
#' calculation of life expectancy from reproductive maturity (i.e., first
#' reproduction), they should account for the possibility that there may be
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
#' # mean life expectancy starting from stage class 2 
#' life_expect_mean(mpm1$matU, start = 2)
#' life_expect_mean(mpm1$matU, start = "small")  # equivalent using named life stages
#' 
#' # mean life expectancy starting from first reproduction
#' rep_stages <- repro_stages(mpm1$matF)
#' n1 <- mature_distrib(mpm1$matU, start = 2, repro_stages = rep_stages)
#' life_expect_mean(mpm1$matU, start = n1)
#'
#'# variance of life expectancy from stage class 1
#' life_expect_var(mpm1$matU, start = 1)
#'
#' @rdname life_expect
#' @export life_expect_mean
life_expect_mean <- function(matU, start = 1L) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU, start_vec = TRUE)
  
  # matrix dimension
  matDim <- nrow(matU)
  
  if (length(start) > 1) {
    start_vec <- start / sum(start)
  } else {
    start_vec <- rep(0.0, matDim)
    if(!is.null(dimnames(matU))) {
      checkMatchingStageNames(matU)
      names(start_vec) <- colnames(matU)
    }
    start_vec[start] <- 1.0
  }
  
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)
  
  if(inherits(N, "try-error")) {
    mean <- NA_real_
  } else {

    mean <- sum(colSums(N) * start_vec)
  }
  
  
  life_expect_mean <-  mean
  
  return(life_expect_mean)
}

#' @rdname life_expect
#' @export life_expect_var
life_expect_var <- function(matU, start = 1L) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU, start_vec = TRUE)
  
  # matrix dimension
  matDim <- nrow(matU)
  
  if (length(start) > 1) {
    start_vec <- start / sum(start)
  } else {
    start_vec <- rep(0.0, matDim)
    if(!is.null(dimnames(matU))) {
      checkMatchingStageNames(matU)
      names(start_vec) <- colnames(matU)
    }
    start_vec[start] <- 1.0
  }
  
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)
  
  if(inherits(N, "try-error")) {
    var <- NA_real_
  } else {
    
    Nvar <- try(sum(2*N^2-N)-colSums(N)*colSums(N))
    var <- sum(Nvar * start_vec)
    
  }
  
  
  life_expect_var <-  var
  
  return(life_expect_var)
}

#This is the old function from v.0.1.0
#Deprecated now, but retaining for backwards
#compatibility.

#' @rdname life_expect
#' @export life_expect
life_expect <- function(matU, start = 1L) {
  .Deprecated("life_expect_mean")
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU, start_vec = TRUE)
  
  # matrix dimension
  matDim <- nrow(matU)
  
  if (length(start) > 1) {
    start_vec <- start / sum(start)
  } else {
    start_vec <- rep(0.0, matDim)
    if(!is.null(dimnames(matU))) {
      checkMatchingStageNames(matU)
      names(start_vec) <- colnames(matU)
    }
    start_vec[start] <- 1.0
  }
  
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)
  
  if(inherits(N, "try-error")) {
    mean <- NA_real_
    var <- NA_real_
  } else {
    
    Nvar <- try(sum(2*N^2-N)-colSums(N)*colSums(N))
    mean <- sum(colSums(N) * start_vec)
    var <- sum(Nvar * start_vec)
    
  }
  
  
  life_expect <- data.frame("mean" = mean,
                            "var" = var)
  
  return(life_expect)
}

