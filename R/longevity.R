#' Calculate longevity from a matrix population model
#'
#' Calculate longevity (the age \emph{x} at which survivorship for a synthetic cohort falls below some 
#' critical proportion) from a matrix population model
#'
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.,
#'   progression, stasis, and retrogression). Optionally with named rows and
#'   columns indicating the corresponding life stage names.
#' @param start The index (or stage name) of the first stage at which the author
#'   considers the beginning of life. Defaults to \code{1}. Alternately, a numeric vector
#'   giving the starting population vector (in which case \code{length(start)}
#'   must match \code{ncol(matU))}. See section \emph{Starting from multiple stages}.
#' @param x_max The maximum age, in units of the MPM projection interval, to
#'   which survivorship will be calculated. Defaults to \code{1000}.
#' @param lx_crit Proportion of initial cohort remaining before all are considered
#' dead (a value between 0 and 1). Defaults to \code{0.01}.
#' 
#' @return Returns longevity, the integer age at which expected survivorship
#'   falls below \code{lx_crit}. If survivorship doesn't reach \code{lx_crit} by
#'   \code{x_max}, returns \code{NA} and prints a warning message.
#' 
#' @section Starting from multiple stages:
#' Rather than specifying argument \code{start} as a single stage class from
#' which all individuals start life, it may sometimes be desirable to allow for
#' multiple starting stage classes. For example, if we want to start our
#' calculation of longevity from reproductive maturity (i.e., first
#' reproduction), we should account for the possibility that there may be
#' multiple stage classes in which an individual could first reproduce.
#' 
#' To specify multiple starting stage classes, specify argument \code{start} as
#' the desired starting population vector, giving the proportion
#' of individuals starting in each stage class (the length of \code{start}
#' should match the number of columns in the relevant MPM).
#' 
#' @note Note that the units of time in returned values are the same as the
#'   (\code{ProjectionInterval}) of the MPM.
#'   
#' @seealso 
#' \code{\link{mature_distrib}} for calculating the proportion of
#' individuals achieving reproductive maturity in each stage class.
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
#'   Morris, W. F. & Doak, D. F. 2003. Quantitative Conservation Biology:
#'   Theory and Practice of Population Viability Analysis. Sinauer Associates,
#'   Sunderland, Massachusetts, USA
#' 
#' @examples
#' data(mpm1)
#' 
#' longevity(mpm1$matU, start = 2)
#' longevity(mpm1$matU, start = "small")  # equivalent using named life stages
#' longevity(mpm1$matU, start = 2, lx_crit = 0.05)
#' 
#' # starting from first reproduction
#' repstages <- repro_stages(mpm1$matF)
#' n1 <- mature_distrib(mpm1$matU, start = 2, repro_stages = repstages)
#' longevity(mpm1$matU, start = n1)
#' 
#' @export longevity
longevity <- function(matU, start = 1L, x_max = 1000, lx_crit = 0.01) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU, start_vec = TRUE)
  if (lx_crit < 0 | lx_crit > 1) {
    stop("lx_crit must be a proportion between 0 and 1", call. = FALSE)
  }
  
  if (length(start) == 1) {
    start_vec <- rep(0.0, nrow(matU))
    if(!is.null(dimnames(matU))) {
      checkMatchingStageNames(matU)
      names(start_vec) <- colnames(matU)
    }
    start_vec[start] <- 1.0
  } else {
    start_vec <- start
  }
  
  lx <- sum(start_vec)
  t <- 0L
  
  while (lx > lx_crit & t < x_max) {
    start_vec <- matU %*% start_vec
    lx <- sum(start_vec)
    t <- t + 1L
  }
  
  if (lx <= lx_crit) {
    longevity <- t
  } else {
    longevity <- NA_integer_
    warning("survivorship did not reach 'lx_crit' by 'x_max': returning NA")
  }
  
  return(longevity)
}
