#' Derive stage-specific vital rates from a matrix population model
#' 
#' @description 
#' Derive a vector of stage-specific vital rates of survival, growth, shrinkage,
#' stasis, dormancy, or reproduction from a matrix population model. These
#' functions include optional arguments for excluding certain stage classes from
#' the calculation (see \emph{Excluding stages}), and defining the set of
#' biologically-possible transitions (see \emph{Possible transitions}).
#' 
#' This decomposition assume that all transition rates are products of a
#' stage-specific survival term (column sums of \code{matU}) and a lower level
#' vital rate that is conditional on survival (growth, shrinkage, stasis,
#' dormancy, or reproduction). Reproductive vital rates that are not conditional
#' on survival (i.e., within a stage class from which there is no survival) are
#' also allowed.
#' 
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix only containing survival-related transitions; 
#'   progression, stasis, and retrogression).
#' @param matR The reproductive component of a matrix population model (i.e., a
#'   square projection matrix only reflecting transitions due to reproduction; either
#'   sexual, clonal, or both).
#' @param posU A logical matrix of the same dimension as \code{matU}, with
#'   elements indicating whether a given \code{matU} transition is possible
#'   (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matU > 0} (see
#'   Details).
#' @param posR A logical matrix of the same dimension as \code{matR}, with
#'   elements indicating whether a given \code{matR} transition is possible
#'   (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matR > 0} (see
#'   Details).
#' @param exclude Integer, character or logical vector indicating stages for 
#'   which transitions both \emph{to} and \emph{from} the stage should be 
#'   excluded from the calculation of vital rates. See section 
#'   \emph{Excluding stages}.
#' @param exclude_row Integer, character or logical vector indicating stages for 
#'   which transitions both \emph{to} and \emph{from} the stage should be 
#'   excluded from the calculation of vital rates. See section 
#'   \emph{Excluding stages}.
#' @param exclude_col Integer, character or logical vector indicating stages for 
#'   which transitions both \emph{to} and \emph{from} the stage should be 
#'   excluded from the calculation of vital rates. See section 
#'   \emph{Excluding stages}.
#' @param dorm_stages Integer or character vector indicating dormant stage 
#'   classes.
#' @param weights_row Vector of stage-specific weights to apply while summing
#'   vital rates across rows within columns (e.g., reproductive value vector).
#' @param surv_only_na If there is only one possible \code{matU} transition in a
#'   given column, should that transition be attributed exclusively to survival?
#'   If \code{TRUE}, the vital rate of growth/stasis/shrinkage in that column
#'   will be coerced to \code{NA}. If \code{FALSE}, dividing the single
#'   transition by the stage-specific survival probability will always yield a
#'   value of \code{1}. Defaults to \code{TRUE}.
#' 
#' @return Vector of vital rates. Vital rates corresponding to impossible
#'   transitions are coerced to \code{NA} (see \emph{Possible transitions}).
#' 
#' @section Possible transitions:
#' A transition rate of \code{0} within a matrix population model may indicate
#' that the transition is not possible in the given life cycle (e.g., tadpoles
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
#' certain vital rates. For instance, a user might not believe that 'growth' to 
#' a dormant stage class really reflects biological growth, in which case the user 
#' could exclude transitions \emph{to} the dormant stage class using the argument
#' \code{exclude_row}. The user may or may not want to ignore 'growth' transitions
#' \emph{from} the dormant stage class, which can be done using
#' \code{exclude_col}. The argument \code{exclude_col} effectively just coerces
#' the respective vital rate to \code{NA}, to prevent it from getting used in
#' subsequent calculations. To exclude transitions both \emph{to and from} a
#' given set of stages, use argument \code{exclude}.
#' 
#' @author Patrick Barks <patrick.barks@@gmail.com>
#'
#' @family vital rates
#' 
#' @examples
#' # create example MPM (stage 4 is dormant)
#' matU <- rbind(c(0.1,   0,   0,   0),
#'               c(0.5, 0.2, 0.1, 0.1),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.4))
#' 
#' matR <- rbind(c(  0,   0.7, 1.1, 0),
#'               c(  0,   0.3, 0.8, 0),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#' 
#' vr_vec_survival(matU, exclude_col = 4)
#' vr_vec_growth(matU, exclude = 4)
#' 
#' # `exclude*` and `*_stages` arguments can accept stage names
#' matU <- name_stages(matU)
#' matR <- name_stages(matR)
#' vr_vec_shrinkage(matU, exclude = 4)
#' vr_vec_stasis(matU, exclude = "stage_4")
#' 
#' vr_vec_dorm_enter(matU, dorm_stages = 4)
#' vr_vec_dorm_exit(name_stages(matU), dorm_stages = "stage_4")
#' 
#' vr_vec_reproduction(matU, matR, exclude_col = "stage_4")
#' 
#' @name vr_vec
NULL


#' @rdname vr_vec
#' @export vr_vec_survival
vr_vec_survival <- function(matU,
                            posU = matU > 0,
                            exclude_col = NULL) {
  
  checkValidMat(matU)
  checkValidStages(matU, exclude_col)
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
                          exclude = NULL,
                          exclude_row = NULL,
                          exclude_col = NULL,
                          surv_only_na = TRUE) {
  
  checkValidMat(matU)
  checkValidStages(matU, exclude)
  checkValidStages(matU, exclude_row)
  checkValidStages(matU, exclude_col)
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU,
                   surv_only_na = surv_only_na)
  
  tri_low <- lower.tri(vmat, diag = FALSE)
  vmat[!tri_low] <- NA_real_
  
  vmat[exclude, ] <- NA_real_
  vmat[exclude_row, ] <- NA_real_
  
  v <- colSums2(vmat)
  
  v[exclude] <- NA_real_
  v[exclude_col] <- NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_shrinkage
vr_vec_shrinkage <- function(matU,
                             posU = matU > 0,
                             exclude = NULL,
                             exclude_row = NULL,
                             exclude_col = NULL,
                             surv_only_na = TRUE) {
  
  checkValidMat(matU)
  checkValidStages(matU, exclude)
  checkValidStages(matU, exclude_row)
  checkValidStages(matU, exclude_col)
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU,
                   surv_only_na = surv_only_na)
  
  tri_upp <- upper.tri(vmat, diag = FALSE)
  vmat[!tri_upp] <- NA_real_
  
  vmat[exclude, ] <- NA_real_
  vmat[exclude_row, ] <- NA_real_
  
  v <- colSums2(vmat)
  
  v[exclude] <- NA_real_
  v[exclude_col] <- NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_stasis
vr_vec_stasis <- function(matU,
                          posU = matU > 0,
                          exclude = NULL,
                          surv_only_na = TRUE) {
  
  checkValidMat(matU)
  checkValidStages(matU, exclude)
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU,
                   surv_only_na = surv_only_na)
  
  v <- diag(vmat)
  v[exclude] <- NA_real_
  
  pos_vital <- diag(posU)
  v[!pos_vital] <- NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_dorm_enter
vr_vec_dorm_enter <- function(matU,
                              posU = matU > 0,
                              dorm_stages) {
  
  checkValidMat(matU)
  checkValidStages(matU, dorm_stages)
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU)
  
  v <- colSums2(vmat[dorm_stages, , drop = FALSE])
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
  
  checkValidMat(matU)
  checkValidStages(matU, dorm_stages)
  
  vmat <- vr_mat_U(matU = matU,
                   posU = posU)
  
  vmat[dorm_stages, ] <- 0
  v <- colSums2(vmat)
  
  # convert dorm_stages to indices, if necessary
  if (is.character(dorm_stages)) {
    dorm_stages <- which(colnames(matU) == dorm_stages)
  }
  
  # possible transitions from dormant to non-dormant
  pos_exit <- matrix(FALSE, nrow = nrow(matU), ncol = nrow(matU))
  pos_exit[-dorm_stages, dorm_stages] <- TRUE
  pos_vital <- apply(posU & pos_exit, 2, any)
  v[!pos_vital] <- NA_real_
  
  return(v)
}


#' @rdname vr_vec
#' @export vr_vec_reproduction
vr_vec_reproduction <- function(matU,
                             matR,
                             posR = matR > 0,
                             exclude_col = NULL,
                             weights_row = NULL) {
  
  checkValidMat(matU)
  checkValidMat(matR)
  checkMatchingStageNames(matU, matR)
  checkValidStages(matU, exclude_col)
  
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

