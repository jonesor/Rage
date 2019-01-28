#' Derive survival-independent vital rates for growth, stasis, shrinkage, and
#' reproduction
#' 
#' @description 
#' Divides columns of a matrix population model by the corresponding
#' stage-specific survival probability, to obtain lower-level vital rates for
#' growth, stasis, shrinkage, and reproduction. Vital rates corresponding to
#' biologically impossible transitions are coerced to \code{NA}.
#'
#' With one exception, these functions assume that transition rates in the
#' matrix population model were in fact calculated as the product of
#' stage-specific survival values and lower-level vital rates of growth, stasis,
#' shrinkage, and reproduction. The one exception is that, if a matrix
#' population model has non-zero reproduction in a stage from which there is no
#' survival, the full reproductive transition will be returned as the vital rate
#' (i.e. because there is no survival component to pull out).
#' 
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @param posU A logical matrix of the same dimension as \code{matU}, with
#'   elements indicating whether a given \code{matU} transition is possible
#'   (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matU > 0} (see
#'   Details).
#' @param posR A logical matrix of the same dimension as \code{matR}, with
#'   elements indicating whether a given \code{matR} transition is possible
#'   (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matR > 0} (see
#'   Details).
#' @param surv_only_na If there is only one possible \code{matU} transition in a
#'   given column, should that transition be attributed exclusively to survival?
#'   If \code{TRUE}, the vital rate of growth/stasis/shrinkage in that column
#'   will be coerced to \code{NA}. If \code{FALSE}, dividing the single
#'   transition by the stage-specific survival probability will always yield a
#'   value of \code{1}. Defaults to \code{TRUE}.
#' 
#' @details 
#' A transition rate of \code{0} within a matrix population model may indicate
#' that the transition is not possible in the given life cycle (e.g. tadpoles
#' never revert to eggs), or that the transition is possible but was estimated
#' to be \code{0} in the relevant population and time period. If vital rates are
#' to be averaged across multiple stage classes, or compared across populations,
#' it may be important to distinguish between these two types of zeros.
#' 
#' By default, the \code{vitals_mat} functions assume that a transition rate of
#' \code{0} indicates an impossible transition, in which case a value of
#' \code{NA} will be returned in the relevant matrix cell. Specifically, the
#' arguments \code{posU} and \code{posR} are specified by the logical
#' expressions \code{(matU > 0)} and \code{(matR > 0)}, respectively. If the
#' matrix population model includes transitions that are possible but estimated
#' to be \code{0}, one should specify the \code{posU} and/or \code{posR}
#' arguments manually.
#' 
#' @return A matrix of vital rates. Vital rates corrsponding to impossible
#'   transitions will be coerced to \code{NA} (see Details).
#' 
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' 
#' @examples
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#' 
#' matR <- rbind(c(  0,   0, 1.1, 1.6),
#'               c(  0,   0, 0.8, 0.4),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#' 
#' # extract vital rates of survival from matU
#' vitals_matU(matU)
#' 
#' # extract vital rates of survival from matR
#' vitals_matR(matU, matR)
#' 
#' @name vitals_mat


#' @rdname vitals_mat
#' @export vitals_matU
vitals_matU <- function(matU, posU = matU > 0, surv_only_na = TRUE) {
  
  checkValidMat(matU)
  sigma <- colSums(matU, na.rm = TRUE)
  sigma[sigma == 0] <- 1 # avoid NaN if transition possible but not observed
  
  vmat <- t(t(matU) / sigma)
  vmat[!posU] <- NA_real_
  
  if (surv_only_na) {
    surv_only <- colSums(posU) == 1
    vmat[,surv_only] <- NA_real_
  }
  
  return(vmat)
}


#' @rdname vitals_mat
#' @export vitals_matR
vitals_matR <- function(matU, matR, posR = matR > 0) {
  
  checkValidMat(matU)
  checkValidMat(matR)
  sigma <- colSums(matU)
  
  sigma[sigma == 0] <- 1 # avoid NaN if no survival
  
  vmat <- t(t(matR) / sigma)
  vmat[!posR] <- NA_real_
  return(vmat)
}

