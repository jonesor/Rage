#' Calculate net reproductive value (R0) from a matrix population model
#'
#' Calculate net reproductive value (R0) from a matrix population model.
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @param start Index of the first stage at which the author considers the
#'   beginning of life. Only used if \code{method = "start"}. Defaults to 1.
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
#' net_repro_rate(matU, matF)
#' 
#' # calculate R0 using the start method
#' net_repro_rate(matU, matF, method = "start", start = 2)
#' 
#' @importFrom popbio lambda
#' @export net_repro_rate
net_repro_rate <- function(matU, matR, start = 1, method = "generation") {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidMat(matR)
  checkValidStartLife(start, matU)
  method <- match.arg(method, c("generation", "start"))
                 
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
                 "start" = R[start, start])
  }
  
  return(R0)
}
