#' Derive stage-specific vital rates from a matrix population model
#' 
#' @description 
#' Derive stage-specific vital rates from a matrix population model. Available
#' vital rates include survival, and a variety of traits that are conditional on
#' survival including growth, shrinkage, stasis, entering and exiting dormancy,
#' and reproduction. Vital rates corresponding to impossible transitions are
#' coerced to \code{NA}.
#' 
#' With one exception, these functions assume that transition rates in the
#' matrix population model were in fact calculated as the product of
#' stage-specific survival values and lower-level vital rates of growth, stasis,
#' shrinkage, and reproduction. The one exception is that, if a matrix
#' population model has non-zero reproduction in a stage from which there is no
#' survival, \code{vr_scalar_fecund} will return the full reproductive transition
#' as the vital rate (i.e. because there is no survival component to pull out).
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
#' @param exclude_row Integer or logical vector indicating stages for which
#'   transitions \emph{to} the stage should be excluded from the calculation of
#'   vital rates. See section 'Excluding stages'.
#' @param exclude_col Integer or logical vector indicating stages for which
#'   transitions \emph{from} the stage should be ignore (coerced to \code{NA}).
#'   See section \emph{Excluding stages}.
#' @param dorm_stages Integer indicator(s) for dormant stage classes.
#' @param weights_row Vector of stage-specific weights to apply while averaging
#'   vital rates across rows within columns. See section \emph{Weighting}.
#' @param weights_col Vector of stage-specific weights to apply while averaging
#'   vital rates across columns. See section \emph{Weighting}.
#' @param surv_only_na If there is only one possible \code{matU} transition in a
#'   given column, should that transition be attributed exclusively to survival?
#'   If \code{TRUE}, the vital rate of growth/stasis/shrinkage in that column
#'   will be coerced to \code{NA}. If \code{FALSE}, dividing the single
#'   transition by the stage-specific survival probability will always yield a
#'   value of \code{1}. Defaults to \code{TRUE}.
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
#' By default, the \code{vr_scalar_} functions assume that a transition rate of
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
#' vr_scalar_survival(matU)
#' vr_scalar_growth(matU)
#' vr_scalar_shrinkage(matU)
#' vr_scalar_stasis(matU)
#' 
#' vr_scalar_dorm_enter(matU, dorm_stages = 4)
#' vr_scalar_dorm_exit(matU, dorm_stages = 4)
#' 
#' vr_scalar_fecund(matU, matF)
#' 
#' @name vr_scalar
NULL


#' @rdname vr_scalar
#' @export vr_scalar_survival
vr_scalar_survival <- function(matU,
                               posU = matU > 0,
                               weights_col = NULL,
                               exclude_col = NULL) {
  
  vr_vec <- vr_vec_survival(matU,
                            posU = posU,
                            exclude_col = exclude_col)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr_scalar
#' @export vr_scalar_growth
vr_scalar_growth <- function(matU,
                             posU = matU > 0,
                             exclude_row = NULL,
                             exclude_col = NULL,
                             weights_col = NULL,
                             surv_only_na = TRUE) {
  
  vr_vec <- vr_vec_growth(matU = matU,
                          posU = posU,
                          exclude_row = exclude_row,
                          exclude_col = exclude_col,
                          surv_only_na = surv_only_na)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr_scalar
#' @export vr_scalar_shrinkage
vr_scalar_shrinkage <- function(matU,
                                posU = matU > 0,
                                weights_col = NULL,
                                exclude_row = NULL,
                                exclude_col = NULL,
                                surv_only_na = TRUE) {
  
  vr_vec <- vr_vec_shrinkage(matU = matU,
                             posU = posU,
                             exclude_row = exclude_row,
                             exclude_col = exclude_col,
                             surv_only_na = surv_only_na)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr_scalar
#' @export vr_scalar_stasis
vr_scalar_stasis <- function(matU,
                             posU = matU > 0,
                             weights_col = NULL,
                             exclude_col = NULL,
                             surv_only_na = TRUE) {
  
  vr_vec <- vr_vec_stasis(matU,
                          posU = posU,
                          exclude_col = exclude_col,
                          surv_only_na = surv_only_na)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr_scalar
#' @export vr_scalar_dorm_enter
vr_scalar_dorm_enter <- function(matU,
                                 posU = matU > 0,
                                 dorm_stages,
                                 weights_col = NULL) {
  
  vr_vec <- vr_vec_dorm_enter(matU = matU,
                              posU = posU,
                              dorm_stages = dorm_stages)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr_scalar
#' @export vr_scalar_dorm_exit
vr_scalar_dorm_exit <- function(matU,
                                posU = matU > 0,
                                dorm_stages,
                                weights_col = NULL) {
  
  vr_vec <- vr_vec_dorm_exit(matU = matU,
                             posU = posU,
                             dorm_stages = dorm_stages)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr_scalar
#' @export vr_scalar_fecund
vr_scalar_fecund <- function(matU,
                             matR,
                             posR = matR > 0,
                             weights_row = NULL,
                             weights_col = NULL) {
  
  vr_vec <- vr_vec_fecund(matU = matU,
                          matR = matR,
                          posR = posR,
                          weights_row = weights_row)
  
  return(column_weight(vr_vec, weights_col))
}



#' @noRd
column_weight <- function(x, weights) {
  if (is.null(weights)) {
    out <- mean(x, na.rm = TRUE)
  } else {
    weights[is.na(x)] <- NA_real_
    weights <- weights / sum(weights, na.rm = TRUE)
    out <- sum(weights * x, na.rm = TRUE)
  }
  return(out)
}
