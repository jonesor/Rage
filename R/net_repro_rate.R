#' Calculate net reproductive rate (R0) from a matrix population model
#'
#' @description
#' Calculate net reproductive rate (R0) from a matrix population model. The net
#' reproduction rate (R0) is the mean number of recruits produced during the
#' mean life expectancy of an individual. See section 5.3.5 of Caswell (2001).
#'
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression). Optionally with named rows and
#'   columns indicating the corresponding life stage names.
#' @param matR The reproductive component of a matrix population model (i.e., a
#'   square projection matrix only reflecting transitions due to reproduction;
#'   either sexual, clonal, or both). If \code{matR} is not provided, it will be
#'   constructed by summing \code{matF} and \code{matC}.
#' @param matF The matrix reflecting sexual reproduction. If provided
#'   without \code{matC}, \code{matC} is assumed to be a zero matrix. If
#'   \code{matR} is provided, this argument is ignored.
#' @param matC The matrix reflecting clonal (asexual) reproduction.
#'   If provided without \code{matF}, \code{matF} is assumed to be a zero
#'   matrix. If \code{matR} is provided, this argument is ignored.
#' @param start Index (or stage name) of the first stage at which the author
#'   considers the beginning of life. Only used if \code{method = "start"}.
#'   Defaults to \code{1}.
#' @param method The method used to calculate net reproductive rate, either
#'   \code{"generation"} or \code{"start"}. Defaults to \code{"generation"}.
#'   See Details.
#'
#' @details
#' The \code{method} argument controls how net reproductive rate is calculated.
#'
#' If \code{method = "generation"}, net reproductive rate is calculated as the
#' per-generation population growth rate (i.e., the dominant eigenvalue of
#' \code{matR \%*\% N}, where \code{N} is the fundamental matrix). See Caswell
#' (2001) Section 5.3.4.
#'
#' If \code{method = "start"}, net reproductive rate is calculated as the
#' expected lifetime production of offspring that start life in stage
#' \code{start}, by an individual also starting life in stage \code{start}
#' (i.e., \code{(matR \%*\% N)[start,start]}).
#'
#' If offspring only arise in stage \code{start}, the two methods give the
#' same result.
#' @return Returns the net reproductive rate. If \code{matU} is singular (often
#'   indicating infinite life expectancy), returns \code{NA}.
#'
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
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
#' net_repro_rate(mpm1$matU, mpm1$matF)
#'
#' # calculate R0 using the start method, specifying either the life stage index
#' # or name
#' net_repro_rate(mpm1$matU, mpm1$matF, method = "start", start = 1)
#' net_repro_rate(mpm1$matU, mpm1$matF, method = "start", start = "seed")
#' 
#' # It is usually better to explicitly name the arguments, rather than relying
#' # on order.
#' net_repro_rate(matU = mpm1$matU, matF = mpm1$matF, 
#' method = "start", start = 1)
#' 
#' net_repro_rate(matU = mpm1$matU, matR = mpm1$matF, 
#' method = "start", start = "seed")
#'
#' @export net_repro_rate
net_repro_rate <- function(matU, matR = NULL, matF = NULL, 
                           matC = NULL, start = 1, method = "generation") {
  
  # Call the helper function to construct matR if not provided
  matR <- process_fertility_inputs(matR, matF, matC)
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidMat(matR)
  checkValidStartLife(start, matU)
  method <- match.arg(method, c("generation", "start"))
  if (!is.numeric(start)) {
    checkMatchingStageNames(M = matU, N = matR)
  }

  # matrix dimensions
  matDim <- nrow(matU)

  # try calculating fundamental matrix (will fail if matrix singular)
  N <- try(solve(diag(matDim) - matU), silent = TRUE)

  # calculate R0
  # first check for singular matU (if singular, R0 = NA)
  if (inherits(N, "try-error") && grepl("singular", N[1], fixed = TRUE)) {
    R0 <- NA_real_
  } else {
    R <- matR %*% N

    R0 <- switch(method,
      "generation" = lambda(R),
      "start" = R[start, start]
    )
  }

  return(R0)
}
