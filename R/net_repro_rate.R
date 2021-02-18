#' Calculate net reproductive value (R0) from a matrix population model
#' 
#' @description 
#' Calculate net reproductive value (R0) from a matrix population model. The net 
#' reproduction value (R0) is the mean number of recruits produced during the 
#' mean life expectancy of an individual. See section 5.3.5 of Caswell (2001).
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression). Optionally with named rows and
#'   columns indicating the corresponding life stage names.
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix only reflecting transitions due to reproduction; either
#'   sexual, clonal, or both). Optionally with named rows and columns indicating
#'   the corresponding life stage names.
#' @param start Index (or stage name) of the first stage at which the author
#'   considers the beginning of life. Only used if \code{method = "start"}.
#'   Defaults to 1.
#' @param method The method used to calculate net reproductive value, either
#'   \code{"generation"} or \code{"start"}. Defaults to \code{"generation"}.
#'   See Details.
#' 
#' @details
#' The \code{method} argument controls how net reproductive rate is calculated.
#' 
#' If \code{method = "generation"}, net reproductive value is calculated as the
#' per-generation population growth rate (i.e. the dominant eigenvalue of
#' \code{matR \%*\% N}, where \code{N} is the fundamental matrix). See Caswell
#' (2001) Section 5.3.4.
#' 
#' If \code{method = "start"}, net reproductive value is calculated as the
#' expected lifetime production of offspring that start life in stage
#' \code{start}, by an individual also starting life in stage \code{start} (i.e.
#' \code{(matR \%*\% N)[start,start]}).
#' 
#' If offspring only arise in stage \code{start}, the two methods give the
#' same result.
#' @return Returns the net reproductive value. If \code{matU} is singular (often
#'   indicating infinite life expectancy), returns \code{NA}.
#'
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
#' 
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' 
#' @examples
#' data(mpm1)
#' 
#' net_repro_rate(mpm1$matU, mpm1$matF)
#' 
#' # calculate R0 using the start method, specifying either the life stage index
#' # or name
#' net_repro_rate(mpm1$matU, mpm1$matF, method = "start", start = 2)
#' net_repro_rate(mpm1$matU, mpm1$matF, method = "start", start = "small")
#' 
#' @importFrom popbio lambda
#' @export net_repro_rate
net_repro_rate <- function(matU, matR, start = 1, method = "generation") {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidMat(matR)
  checkValidStartLife(start, matU)
  method <- match.arg(method, c("generation", "start"))
  if (!is.numeric(start)){
    checkMatchingStageNames(M = matU, N = matR)
  }
  
  # matrix dimensions
  matDim <- nrow(matU)
  
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)
  
  # calculate R0
  # first check for singular matU (if singular, R0 = NA)
  if (("try-error" %in% class(N)) && grepl('singular', N[1])) {
    R0 <- NA_real_
  } else {
    R <- matR %*% N
    
    R0 <- switch(method,
                 "generation" = popbio::lambda(R),
                 "start" = R[start, start])
  }
  
  return(R0)
}
