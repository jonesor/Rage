#' Calculate net reproductive value from a matrix population model
#'
#' Calculate net reproductive value from a matrix population model.
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @param startLife Index of the first stage at which the author considers the
#'   beginning of life. Only used if \code{method = "startLife"}. Defaults to 1.
#' @param method The method used to calculate net reproductive value, either
#'   \code{"generation"} or \code{"startLife"}. Defaults to \code{"generation"}.
#'   See Details.
#' @details
#' The \code{method} argument controls how net reproductive rate is calculated.
#' 
#' If \code{method = "generation"}, net reproductive value is calculated as the
#' per-generation population growth rate (i.e. the dominant eigenvalue of
#' \code{matR \%*\% N}, where \code{N} is the fundamental matrix). See Caswell
#' (2001) Section 5.3.4.
#' 
#' If \code{method = "startLife"}, net reproductive value is calculated as the
#' expected lifetime production of offspring that start life in stage
#' \code{startLife}, by an individual also starting life in stage
#' \code{startLife} (i.e. \code{(matR \%*\% N)[startLife,startLife]}).
#' 
#' If offspring only arise in stage \code{startLife}, the two methods give the
#' same result.
#' @return Returns the net reproductive value. If \code{matU} is singular (often
#'   indicating infinite life expectancy), returns \code{NA}.
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' @examples
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#' 
#' matF <- rbind(c(  0,   0, 1.1, 1.6),
#'               c(  0,   0, 0.8, 0.4),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#' 
#' R0(matU, matF)
#' 
#' # calculate R0 using the startLife method
#' R0(matU, matF, method = "startLife", startLife = 2)
#' 
#' @importFrom popbio lambda
#' @export R0
R0 <- function(matU, matR, startLife = 1, method = "generation") {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidMat(matR)
  checkValidStartLife(startLife, matU)
  method <- match.arg(method, c("generation", "startLife"))
                 
  # matrix dimensions
  matDim <- nrow(matU)
  
  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)
  
  # calculate R0
  # first check for singular matU (if singular, R0 = NA)
  if (class(N) == 'try-error' && grepl('singular', N[1])) {
    R0 <- NA_real_
  } else {
    R <- matR %*% N
    
    R0 <- switch(method,
                 "generation" = popbio::lambda(R),
                 "startLife" = R[startLife, startLife])
  }
  
  return(R0)
}
