#' Derive stage-specific vital rates from a matrix population model
#' 
#' Derive stage-specific vital rates from a matrix population model. Available
#' vital rates include survival, and a variety of traits that are conditional on
#' survival including growth, shrinkage, stasis, entering and exiting dormancy,
#' and reproduction. Vital rates corresponding to impossible transitions are
#' coerced to \code{NA} (see Details).
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
#' @param stage_exclude Integer indicator(s) for stages that should be excluded
#'   from the calculation of vital rates. Transitions both to and from the
#'   indicated stage will be set to \code{NA}.
#' @param dorm_stages Integer indicator(s) for dormant stage classes.
#' @param weights Vector of stage-specific weights to apply while averaging
#'   vital rates.
#' @param surv_only_na If there is only one possible \code{matU} transition in a
#'   given column, should that transition be attributed exclusively to survival?
#'   If \code{TRUE}, the vital rate of growth/stasis/shrinkage in that column
#'   will be coerced to \code{NA}. If \code{FALSE}, dividing the single
#'   transition by the stage-specific survival probability will always yield a
#'   value of \code{1}. Defaults to \code{TRUE}.
#' @param surv_zero_na If there are reproductive transitions from a stage in
#'   which the survival probability is \code{0}, should the relevant fecundity
#'   vital rates be coerced to \code{NA} (\code{TRUE}), or should the full
#'   fecundity transitions be returned as vital rates (\code{FALSE}). Defaults
#'   to \code{FALSE}.
#' 
#' @return Vector of vital rates. Vital rates corrsponding to impossible
#'   transitions will be coerced to \code{NA} (see Details).
#' 
#' @details 
#' A transition rate of \code{0} within a matrix population model may indicate
#' that the transition is not possible in the given life cycle (e.g. tadpoles
#' never revert to eggs), or that the transition rate is possible but was
#' estimated to be \code{0} in the relevant population and time period. If vital
#' rates are to be averaged across multiple stage classes, or compared across
#' populations, it may be important to distinguish between these two types of
#' zeros.
#' 
#' By default, the \code{vitals_} functions assume that a transition rate of
#' \code{0} indicates an impossible transition, in which case a value of
#' \code{NA} will be used in relevant calculations. Specifically, the arguments
#' \code{posU} and \code{posR} are specified by the logical expressions
#' \code{(matU > 0)} and \code{(matR > 0)}, respectively. If the matrix
#' population model includes transitions that are estimated to be \code{0} but
#' still in fact possible, one should specify the \code{posU} and/or \code{posR}
#' arguments manually.
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
#' matF <- rbind(c(  0,   0, 1.1, 1.6),
#'               c(  0,   0, 0.8, 0.4),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#' 
#' vitals_survival(matU)
#' vitals_growth(matU)
#' vitals_shrinkage(matU)
#' vitals_stasis(matU)
#' 
#' vitals_dorm_enter(matU, dorm_stages = 4)
#' vitals_dorm_exit(matU, dorm_stages = 4)
#' 
#' vitals_fecund(matU, matF)
#' 
#' @name vitals
NULL


#' @rdname vitals
#' @export vitals_survival
vitals_survival <- function(matU, posU = matU > 0, stage_exclude = NULL) {
  
  checkValidMat(matU)
  pos_vital <- apply(posU, 2, any)
  
  v <- colSums(matU)
  v[!pos_vital] <- NA_real_
  v[stage_exclude] <- NA_real_
  return(v)
}


#' @rdname vitals
#' @export vitals_growth
vitals_growth <- function(matU, posU = matU > 0, stage_exclude = NULL,
                          surv_only_na = TRUE) {
  
  vmat <- vitals_matU(matU, posU = posU, surv_only_na = surv_only_na)
  tri_low <- lower.tri(vmat, diag = FALSE)
  vmat[!tri_low] <- NA_real_
  
  if (!is.null(stage_exclude)) vmat[stage_exclude, ] <- NA_real_
  
  v <- colSums2(vmat)
  v[stage_exclude] <- NA_real_
  return(v)
}


#' @rdname vitals
#' @export vitals_shrinkage
vitals_shrinkage <- function(matU, posU = matU > 0, stage_exclude = NULL,
                             surv_only_na = TRUE) {
  
  vmat <- vitals_matU(matU, posU = posU, surv_only_na = surv_only_na)
  tri_upp <- upper.tri(vmat, diag = FALSE)
  vmat[!tri_upp] <- NA_real_
  
  if (!is.null(stage_exclude)) vmat[stage_exclude, ] <- NA_real_
  
  v <- colSums2(vmat)
  v[stage_exclude] <- NA_real_
  return(v)
}


#' @rdname vitals
#' @export vitals_stasis
vitals_stasis <- function(matU, posU = matU > 0, stage_exclude = NULL,
                          surv_only_na = TRUE) {
  
  vmat <- vitals_matU(matU, posU = posU, surv_only_na = surv_only_na)
  pos_vital <- diag(posU)
  
  v <- diag(vmat)
  v[stage_exclude] <- NA_real_  # set excluded stages to NA
  v[!pos_vital] <- NA_real_     # set stages with no possible stasis to NA
  return(v)
}


#' @rdname vitals
#' @export vitals_dorm_enter
vitals_dorm_enter <- function(matU, posU = matU > 0, dorm_stages) {
  
  vmat <- vitals_matU(matU, posU = posU)
  pos_vital <- apply(posU[dorm_stages, , drop = FALSE], 2, any)
  
  v <- colSums(vmat[dorm_stages, , drop = FALSE], na.rm = TRUE)
  v[dorm_stages] <- NA_real_
  v[!pos_vital] <-  NA_real_
  return(v)
}


#' @rdname vitals
#' @export vitals_dorm_exit
vitals_dorm_exit <- function(matU, posU = matU > 0, dorm_stages) {
  
  vmat <- vitals_matU(matU, posU = posU)
  
  # possible transitions from dormant to non-dormant
  pos_exit <- matrix(FALSE, nrow = nrow(matU), ncol = nrow(matU))
  pos_exit[-dorm_stages, dorm_stages] <- TRUE
  pos_vital <- apply(posU & pos_exit, 2, any)
  
  vmat[dorm_stages, ] <- 0
  v <- colSums(vmat)
  v[!pos_vital] <- NA_real_
  return(v)
}


#' @rdname vitals
#' @export vitals_fecund
vitals_fecund <- function(matU, matR, posR = matR > 0, weights = NULL,
                          surv_zero_na = FALSE) {
  
  vmat <- vitals_matR(matU, matR, posR = posR, surv_zero_na = surv_zero_na)
  pos_vital <- apply(posR, 2, any)
  
  if (is.null(weights)) weights <- rep(1, nrow(matU))
  v <- colSums2(vmat * weights)
  return(v)
}

