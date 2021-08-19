#' Calculate age-specific traits from a matrix population model
#'
#' These functions use age-from-stage decomposition methods to calculate
#' age-specific survivorship (\code{lx}), survival probability (\code{px}), mortality hazard
#' (\code{hx}), or reproduction (\code{mx}) from a matrix population model (MPM). A detailed
#' description of these methods can be found in sections 5.3.1 and 5.3.2 of
#' Caswell (2001).
#' 
#' @param matU The survival component of a MPM (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.,
#'   progression, stasis, and retrogression). Optionally with named rows and
#'   columns indicating the corresponding life stage names.
#' @param matR The reproductive component of a MPM (i.e., a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both). Optionally with named rows and columns indicating
#'   the corresponding life stage names.
#' @param start The index (or stage name) of the first stage at which the author
#'   considers the beginning of life. Defaults to 1. Alternately, a numeric vector
#'   giving the starting population vector (in which case \code{length(start)}
#'   must match \code{ncol(matU))}. See section \emph{Starting from multiple stages}.
#' @param xmax Maximum age to which age-specific traits will be calculated
#'   (defaults to \code{100000}).
#' @param lx_crit Minimum value of \code{lx} to which age-specific traits will be
#'   calculated (defaults to \code{0.0001}).
#' @param tol To account for floating point errors that occasionally lead to
#'   values of \code{lx} slightly greater than 1, values of \code{lx} within the open interval
#'   (\code{1}, \code{1 + tol}) are coerced to 1. Defaults to \code{0.0001}. To
#'   prevent coercion, set \code{tol} to \code{0}.
#' 
#' @return A vector
#' 
#' @note Note that the units of time for the returned vectors (i.e., \code{x}) are
#'   the same as the projection interval (\code{ProjectionInterval}) of the MPM.
#' 
#' 
#' @section Starting from multiple stages:
#' Rather than specifying argument \code{start} as a single stage class from
#' which all individuals start life, it may sometimes be desirable to allow for
#' multiple starting stage classes. For example, if users want to start their
#' calculation of age-specific traits from reproductive maturity (i.e., first
#' reproduction), they should account for the possibility that there may be
#' multiple stage classes in which an individual could first reproduce.
#' 
#' To specify multiple starting stage classes, users should specify argument
#' \code{start} as the desired starting population vector (\strong{n1}), giving
#' the proportion of individuals starting in each stage class (the length of
#' \code{start} should match the number of columns in the relevant MPM).
#' 
#' See function \code{\link{mature_distrib}} for calculating the proportion of
#' individuals achieving reproductive maturity in each stage class.
#' 
#' @note The output vector is calculated recursively until the age class (\code{x})
#'   reaches \code{xmax} or survivorship (\code{lx}) falls below \code{lx_crit},
#'   whichever comes first. To force calculation to \code{xmax}, set
#'   \code{lx_crit} to \code{0}. Conversely, to force calculation to
#'   \code{lx_crit}, set \code{xmax} to \code{Inf}.
#'
#' @note Note that the units of time in returned values (i.e., \code{x}) are the
#'   same as the projection interval (`ProjectionInterval`) of the MPM.
#'   
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @family life tables
#' 
#' @seealso \code{\link{lifetable_convert}}
#' 
#' @references 
#' Caswell, H. 2001. Matrix Population Models: Construction, Analysis, and
#' Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968
#' 
#' Jones O. R. 2021. Life tables: Construction and interpretation In:
#' Demographic Methods Across the Tree of Life. Edited by Salguero-Gomez R &
#' Gamelon M. Oxford University Press. Oxford, UK. ISBN: 9780198838609
#' 
#' Preston, S., Heuveline, P., & Guillot, M. 2000. Demography: Measuring and
#' Modeling Population Processes. Wiley. ISBN: 9781557864512
#' 
#' @examples
#' data(mpm1)
#' 
#' # age-specific survivorship
#' mpm_to_lx(mpm1$matU)
#' mpm_to_lx(mpm1$matU, start = 2)       # starting from stage 2
#' mpm_to_lx(mpm1$matU, start = "small") # equivalent using named life stages
#' mpm_to_lx(mpm1$matU, xmax = 10)       # to a maximum age of 10
#' mpm_to_lx(mpm1$matU, lx_crit = 0.05)  # to a minimum lx of 0.05
#' 
#' # age-specific survival probability
#' mpm_to_px(mpm1$matU)
#' 
#' # age-specific mortality hazard
#' mpm_to_hx(mpm1$matU)
#' 
#' # age-specific fecundity
#' mpm_to_mx(mpm1$matU, mpm1$matF)
#' 
#' 
#' ### starting from first reproduction
#' repstages <- repro_stages(mpm1$matF)
#' n1 <- mature_distrib(mpm1$matU, start = 2, repro_stages = repstages)
#' 
#' mpm_to_lx(mpm1$matU, start = n1)
#' mpm_to_px(mpm1$matU, start = n1)
#' mpm_to_hx(mpm1$matU, start = n1)
#' mpm_to_mx(mpm1$matU, mpm1$matF, start = n1)
#' 
#' @name age_from_stage
NULL


#' @rdname age_from_stage
#' @export mpm_to_mx
mpm_to_mx <- function(matU, matR, start = 1L, xmax = 1e5, lx_crit = 1e-4,
                      tol = 1e-4) {
  
  # validate arguments (leave rest to mpm_to_lx)
  checkValidMat(matR)
  
  N <- length(mpm_to_lx(matU, start, xmax, lx_crit, tol))
  
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
  
  n <- start_vec
  mx <- numeric(N)
  
  for (i in 1:N) {
    n <- n / sum(n)
    phi <- matR %*% n
    mx[i] <- sum(phi)
    n <- matU %*% n
  }
  
  mx[is.nan(mx)] <- 0
  return(mx)
}


#' @rdname age_from_stage
#' @export mpm_to_lx
mpm_to_lx <- function(matU, start = 1L, xmax = 1e5, lx_crit = 1e-4,
                      tol = 1e-4) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  checkValidStartLife(start, matU, start_vec = TRUE)
  
  t <- 0L
  lx <- 1.0
  lx_vec <- lx
  
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
  
  n <- start_vec
  
  while (lx > lx_crit & t < xmax) {
    n <- matU %*% n
    lx <- sum(n)
    t <- t + 1L
    lx_vec[t] <- lx
  }
  
  # fix floating point errors resulting in lx > 1.0
  if (tol > 0) {
    lx_vec[lx_vec > 1.0 & lx_vec < (1.0 + tol)] <- 1.0
  }
  
  return(c(1.0, lx_vec))
}


#' @rdname age_from_stage
#' @export mpm_to_px
mpm_to_px <- function(matU, start = 1L, xmax = 1e5, lx_crit = 1e-4,
                      tol = 1e-4) {
  # leave argument validation to mpm_to_lx
  lx <- mpm_to_lx(matU, start, xmax, lx_crit, tol)
  return(lx_to_px(lx))
}


#' @rdname age_from_stage
#' @export mpm_to_hx
mpm_to_hx <- function(matU, start = 1L, xmax = 1e5, lx_crit = 1e-4,
                      tol = 1e-4) {
  # leave argument validation to mpm_to_lx
  lx <- mpm_to_lx(matU, start, xmax, lx_crit, tol)
  return(lx_to_hx(lx))
}
