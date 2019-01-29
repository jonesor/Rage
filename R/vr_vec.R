#' Derive stage-specific vital rates from a matrix population model
#' 
#' @description 
#' Derive a vector of stage-specific vital rates of survival, growth, shrinkage,
#' stasis, dormancy, or reproduction from a matrix population model. These
#' functions include optional arguments for excluding certain stage classes from
#' the calculation (see \emph{Excluding stages}), and defining the set of
#' biologically-possible transitions (see \emph{Possible transitions}).
#' 
#' These decompositions assume that all transition rates are products of a
#' stage-specific survival term (column sums of \code{matU}) and a lower level
#' vital rate that is conditional on survival (growth, shrinkage, stasis,
#' dormancy, or reproduction). Reproductive vital rates that are not conditional
#' on survival (i.e. within a stage class from which there is no survival) are
#' also allowed.
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
#'   vital rates. See section \emph{Excluding stages}.
#' @param exclude_col Integer or logical vector indicating stages for which
#'   transitions \emph{from} the stage should be ignore (coerced to \code{NA}).
#'   See section \emph{Excluding stages}.
#' @param dorm_stages Integer indicator(s) for dormant stage classes.
#' @param weights_row Vector of stage-specific weights to apply while summing
#'   vital rates across rows within columns (e.g. reproductive value vector).
#' @param surv_only_na If there is only one possible \code{matU} transition in a
#'   given column, should that transition be attributed exclusively to survival?
#'   If \code{TRUE}, the vital rate of growth/stasis/shrinkage in that column
#'   will be coerced to \code{NA}. If \code{FALSE}, dividing the single
#'   transition by the stage-specific survival probability will always yield a
#'   value of \code{1}. Defaults to \code{TRUE}.
#' 
#' @return Vector of vital rates. Vital rates corrsponding to impossible
#'   transitions are coerced to \code{NA} (see \emph{Possible transitions}).
#' 
#' @section Possible transitions:
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
#' @section Excluding stages:
#' It may be desirable to exclude one or more stages from the calculation of
#' certain vital rates. For instance, we might not believe that 'growth' to a
#' dormant stage class really reflects biological growth, in which case we could
#' exclude transitions \emph{to} the dormant stage class using the argument
#' \code{exclude_row}. We may or may not want to ignore 'growth' transitions
#' \emph{from} the dormant stage class, which can be done using
#' \code{exclude_col}. The argument \code{exclude_col} effectively just coerces
#' the respective vital rate to \code{NA}, to prevent it from getting used in
#' subsequent calculations.
#' 
#' @author Patrick Barks <patrick.barks@@gmail.com>
#'   
#' @examples
#' # create example MPM (stage 4 is dormant)
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1, 0.1),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.4))
#' 
#' matF <- rbind(c(  0,   0.7, 1.1, 0),
#'               c(  0,   0.3, 0.8, 0),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#' 
#' vr_vec_survival(matU, exclude_col = 4)
#' vr_vec_growth(matU, exclude_col = 4)
#' vr_vec_shrinkage(matU, exclude_col = 4)
#' vr_vec_stasis(matU, exclude_col = 4)
#' 
#' vr_vec_dorm_enter(matU, dorm_stages = 4)
#' vr_vec_dorm_exit(matU, dorm_stages = 4)
#' 
#' vr_vec_fecundity(matU, matF, exclude_col = 4)
#' 
#' @name vr_vec
NULL


#' @rdname vr_vec
#' @export vr_vec_survival
vr_vec_survival <- function(matU,
                            posU = matU > 0,
                            exclude_col = NULL) {
  
  checkValidMat(matU)
  pos_vital <- apply(posU, 2, any)
  
  v <- colSums(matU)
  v[!pos_vital] <- NA_real_
  v[exclude_col] <- NA_real_
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_growth
vr_vec_growth <- function(matU,
                          posU = matU > 0,
                          exclude_row = NULL,
                          exclude_col = NULL,
                          surv_only_na = TRUE) {
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU,
                   surv_only_na = surv_only_na)
  
  tri_low <- lower.tri(vmat, diag = FALSE)
  vmat[!tri_low] <- NA_real_
  
  vmat[exclude_row, ] <- NA_real_
  
  v <- colSums2(vmat)
  v[exclude_col] <- NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_shrinkage
vr_vec_shrinkage <- function(matU,
                             posU = matU > 0,
                             exclude_row = NULL,
                             exclude_col = NULL,
                             surv_only_na = TRUE) {
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU,
                   surv_only_na = surv_only_na)
  
  tri_upp <- upper.tri(vmat, diag = FALSE)
  vmat[!tri_upp] <- NA_real_
  
  vmat[exclude_row, ] <- NA_real_
  
  v <- colSums2(vmat)
  v[exclude_col] <- NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_stasis
vr_vec_stasis <- function(matU,
                          posU = matU > 0,
                          exclude_col = NULL,
                          surv_only_na = TRUE) {
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU,
                   surv_only_na = surv_only_na)
  
  v <- diag(vmat)
  v[exclude_col] <- NA_real_
  
  pos_vital <- diag(posU)
  v[!pos_vital] <- NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_dorm_enter
vr_vec_dorm_enter <- function(matU,
                              posU = matU > 0,
                              dorm_stages) {
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU)
  
  v <- colSums(vmat[dorm_stages, , drop = FALSE], na.rm = TRUE)
  v[dorm_stages] <- NA_real_
  
  pos_vital <- apply(posU[dorm_stages, , drop = FALSE], 2, any)
  v[!pos_vital] <-  NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_dorm_exit
vr_vec_dorm_exit <- function(matU,
                             posU = matU > 0,
                             dorm_stages) {
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU)
  
  vmat[dorm_stages, ] <- 0
  v <- colSums(vmat)
  
  # possible transitions from dormant to non-dormant
  pos_exit <- matrix(FALSE, nrow = nrow(matU), ncol = nrow(matU))
  pos_exit[-dorm_stages, dorm_stages] <- TRUE
  pos_vital <- apply(posU & pos_exit, 2, any)
  v[!pos_vital] <- NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_fecundity
vr_vec_fecundity <- function(matU,
                             matR,
                             posR = matR > 0,
                             exclude_col = NULL,
                             weights_row = NULL) {
  
  vmat <- vr_mat_R(matU = matU,
                   matR = matR,
                   posR = posR)
  
  if (is.null(weights_row)) weights_row <- rep(1, nrow(matU))
  v <- colSums2(vmat * weights_row)
  
  v[exclude_col] <- NA_real_
  
  pos_vital <- apply(posR, 2, any)
  v[!pos_vital] <- NA_real_
  
  return(v)
}

