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
#'   If set to `NULL` the function returns mean life expectancy from each of the
#'   stages of the MPM.
#' @param mixdist A vector with a length equal to the dimension of the MPM
#'   defining how the function should average the output over the. possible
#'   starting states. See section \emph{Starting from multiple stages}. If this
#'   argument is used, `start` must be set to `NULL`.
#'
#' @section Starting from multiple stages:
#' Sometimes, it is necessary to calculate life expectancy considering multiple
#' starting stage classes instead of just a single stage from which all
#' individuals begin their lives. This scenario arises when there are several
#' possible stages at which an individual can start a particular life event,
#' such as reproductive maturity.
#' To handle such cases, the function provides support for multiple starting
#' stage classes. When calculating life expectancy in this context, the outputs
#' should be averaged using weights determined by the distribution of
#' individuals across these stages. To achieve this, the `start` argument should
#' be set to `NULL`, indicating that the starting stage is not specified, and
#' the `mixdist` argument should be utilized.
#' In the context described, The `mixdist` argument expects a vector that
#' represents the proportion of individuals with their first reproduction in
#' each stage of the MPM. By providing this distribution, the function
#' calculates the mean lifespan by appropriately weighting the life expectancies
#' corresponding to each starting stage.
#' For a practical example that demonstrates this usage, please refer to the
#' code example below.
#'
#' @return Returns life expectancy in the units of the projection interval
#'   (`ProjectionInterval`) of the MPM. If \code{matU} is singular (often
#'   indicating infinite life expectancy), returns \code{NA}.
#'
#' @author Christine M. Hernández <cmh352@cornell.edu>
#' @author Owen R. Jones <jones@biology.sdu.dk>
#'
#' @family life history traits
#'
#' @references Caswell, H. 2001. Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#'
#' @examples
#' data(mpm1)
#'
#' # mean life expectancy starting from stage class 2
#' life_expect_mean(mpm1$matU, start = 2)
#'
#' # equivalent using named life stages
#' life_expect_mean(mpm1$matU, start = "small")
#'
#' # mean life expectancies starting from each of the stages
#' life_expect_mean(mpm1$matU, start = NULL)
#'
#' # mean life expectancy starting from first reproduction, where this varies 
#' across individuals
#' rep_stages <- repro_stages(mpm1$matF)
#' (n1 <- mature_distrib(mpm1$matU, start = 2, repro_stages = rep_stages))
#' life_expect_mean(mpm1$matU, mixdist = n1, start = NULL)
#'
#' # variance of life expectancy from stage class 1
#' life_expect_var(mpm1$matU, start = 1)
#'
#' # variance of life expectancy from stage class 1
#' life_expect_var(mpm1$matU, start = "seed")
#'
#' # variance of life expectancy from each stage class
#' life_expect_var(mpm1$matU, start = NULL)
#'
#' # variance of life expectancies with a set mixing distribution
#' life_expect_var(mpm1$matU, mixdist = c(0.0,0.1,0.3,0.1,0.5), start = NULL)
#' 
#' # setting mixdist to ignore all but one stage should produce the same result
#' as setting the start argument to that stage
#' life_expect_mean(mpm1$matU, start = 3)
#' life_expect_mean(mpm1$matU, mixdist = c(0,0,1,0,0), start = NULL)
#'
#' @rdname life_expect
#' @export life_expect_mean
life_expect_mean <- function(matU, mixdist = NULL, start = 1L) {
  # validate arguments
  # You cannot use both mixdist and start.
  if (!is.null(mixdist) && !is.null(start)) {
    stop("You cannot apply the mixing distribution and also specify a starting
         state. The mixing distribution defines how you want the function to
         average over all possible starting states.")
  }
  
  # check that the MPM is valid
  checkValidMat(matU, warn_surv_issue = TRUE)

  # check that, if it is not NULL, start is valid (i.e. that is it an integer or
  # stage name that matches the MPM)
  if (!is.null(start)) {
    checkValidStartLife(start, matU, start_vec = TRUE)
  }

  # if start is a character string (e.g. a stage name) turn it into a numeric.
  if (inherits(start, "character")) {
    startNumeric <- match(start, colnames(matU))
  } else {
    startNumeric <- start
  }

  matDim <- dim(matU)[1]

  ## Calculate Ex(R | current state)
  # This is an alternative to exactLTRE::fundamental_matrix, to avoid dependency
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)

  # If the calculation of fundamental matrix produces an error, then make the
  # output NA. Otherwise, calculate the life expectancy.
  if (inherits(N, "try-error")) {
    expLCond_z <- NA_real_
  } else {
    expLCond_z <- rep(1, matDim) %*% N
  }


  # If mixdist is NOT null
  if (!is.null(mixdist)) {
    expL <- expLCond_z %*% mixdist
    return(expL)
  }

  # If start is NOT null
  if (!is.null(start)) {
    return(expLCond_z[startNumeric])
  }

  if (is.null(start)) {
    return(expLCond_z)
  }
}

# Calculate the variance in lifespan:
# note: this calculates the variance in the number of time steps!
#' @rdname life_expect
#' @export life_expect_var
life_expect_var <- function(matU, mixdist = NULL, start = 1L) {
  # validate arguments
  # You cannot use both mixdist and start.
  if (!is.null(mixdist) && !is.null(start)) {
    stop("You cannot apply the mixing distribution and also specify a starting
         state. The mixing distribution defines how you want the function to
         average over all possible starting states.")
  }
  
  # check that the MPM is valid
  checkValidMat(matU, warn_surv_issue = TRUE)
  
  # check that, if it is not NULL, start is valid (i.e. that is it an integer or
  # stage name that matches the MPM)
  if (!is.null(start)) {
    checkValidStartLife(start, matU, start_vec = TRUE)
  }
  
  # if start is a character string (e.g. a stage name) turn it into a numeric.
  if (inherits(start, "character")) {
    startNumeric <- match(start, colnames(matU))
  } else {
    startNumeric <- start
  }
  
  matDim <- dim(matU)[1]

  # calculate the fundamental matrix
  # This is an alternative to exactLTRE::fundamental_matrix, to avoid dependency
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)

  if (inherits(N, "try-error")) {
    mean <- NA_real_
  } else {
    expLCond_z <- life_expect_mean(matU, mixdist = NULL, start = NULL)
  }
  ## Calculate Ex(R | current state)

  ## Var(L | current state) using eqn. 5.12 from Hal Caswell (2001)
  eT <- matrix(data = 1, ncol = matDim, nrow = 1) # column vector of 1's
  varLCond_z <- eT %*% (2 * N %*% N - N) - (expLCond_z)^2

  if (is.null(mixdist)) {
    outputVar <- varLCond_z
  } else {
    # variance in LRO due to differences along trajectories:
    varL_within <- varLCond_z %*% mixdist
    # variance in LRO due to differences among starting states:
    varL_between <- t(mixdist) %*% t(expLCond_z^2) -
      (t(mixdist) %*% t(expLCond_z))^2
    # total variance in lifespan, given the mixing distribution:
    outputVar <- varL_within + varL_between
  }
  if(!is.null(start)){
    return(outputVar[startNumeric])
  } 
  if(is.null(start)){
    return(outputVar)
  }
}
