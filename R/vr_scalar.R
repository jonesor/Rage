#' Derive mean vital rates from a matrix population model
#' 
#' @description 
#' Derive mean vital rates of survival, growth, shrinkage, stasis, dormancy, or
#' reproduction from a matrix population model, by averaging across stage
#' classes. These functions include optional arguments for custom weighting of
#' different stage classes (see \emph{Weighting stages}), excluding certain
#' stage classes from the calculation (see \emph{Excluding stages}), and
#' defining the set of biologically-possible transitions (see \emph{Possible
#' transitions}).
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
#'   \emph{Possible transitions}).
#' @param posR A logical matrix of the same dimension as \code{matR}, with
#'   elements indicating whether a given \code{matR} transition is possible
#'   (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{matR > 0} (see
#'   \emph{Possible transitions}).
#' @param exclude_row Integer or logical vector indicating stages for which
#'   transitions \emph{to} the stage should be excluded from the calculation of
#'   vital rates. See section \emph{Excluding stages}.
#' @param exclude_col Integer or logical vector indicating stages for which
#'   transitions \emph{from} the stage should be ignore (coerced to \code{NA}).
#'   See section \emph{Excluding stages}.
#' @param dorm_stages Integer indicator(s) for dormant stage classes.
#' @param weights_row Vector of stage-specific weights to apply while summing
#'   vital rates across rows within columns. See section \emph{Weighting
#'   stages}.
#' @param weights_col Vector of stage-specific weights to apply while averaging
#'   vital rates across columns. See section \emph{Weighting stages}.
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
#' By default, the \code{vr_} functions assume that a transition rate of
#' \code{0} indicates an impossible transition, in which case a value of
#' \code{NA} will be used in relevant calculations. Specifically, the arguments
#' \code{posU} and \code{posR} are specified by the logical expressions
#' \code{(matU > 0)} and \code{(matR > 0)}, respectively. If the matrix
#' population model includes transitions that are estimated to be \code{0} but
#' still in fact possible, one should specify the \code{posU} and/or \code{posR}
#' arguments manually.
#' 
#' @section Weighting stages:
#' In averaging vital rates across stages, it may be desirable to weight stage
#' classes differently (e.g. based on reproductive values, or stable
#' distributions). Weights are generally applied when averaging across columns,
#' i.e., across transitions \emph{from} a set of stage classes (e.g. averaging
#' stage-specific survival probabilities across multiple stages). All \code{vr_}
#' functions therefore include an optional argument \code{weights_from}.
#'
#' In principle, particularly for vital rates of reproduction, we could also
#' apply weights when summing across rows within columns, i.e., across
#' reproductive transitions \emph{to} a set of stage classes (e.g. summing the
#' production of different types of offspring, such as seeds vs. seedlings). For
#' the function \code{vr_fecundity}, we therefore also include an optional
#' argument \code{weights_to}.
#' 
#' If supplied, \code{weights_from} will automatically be scaled to sum to 1
#' over the set of possible transitions, whereas \code{weights_to} will not be
#' rescaled because we wish to enable the use of reproductive values here, which
#' do not naturally sum to 1.
#' 
#' @section Excluding stages:
#' It may be desirable to exclude one or more stages from the calculation of
#' certain vital rates. For instance, we might not believe that 'growth' to a
#' dormant stage class really reflects biological growth, in which case we could
#' exclude transitions \emph{to} the dormant stage class using the argument
#' \code{exclude_row}. We may or may not want to ignore 'growth' transitions
#' \emph{from} the dormant stage class, which can be done using
#' \code{exclude_col}.
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
#' vr_survival(matU, exclude_col = 4)
#' vr_growth(matU, exclude_row = 4, exclude_col = 4)
#' vr_shrinkage(matU, exclude_row = 4, exclude_col = 4)
#' vr_stasis(matU, exclude_col = 4)
#' 
#' vr_dorm_enter(matU, dorm_stages = 4)
#' vr_dorm_exit(matU, dorm_stages = 4)
#' 
#' vr_fecundity(matU, matF, exclude_col = 4)
#' 
#' @name vr
NULL


#' @rdname vr
#' @export vr_survival
vr_survival <- function(matU,
                        posU = matU > 0,
                        exclude_col = NULL,
                        weights_col = NULL) {
  
  vr_vec <- vr_vec_survival(matU = matU,
                            posU = posU,
                            exclude_col = exclude_col)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr
#' @export vr_growth
vr_growth <- function(matU,
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


#' @rdname vr
#' @export vr_shrinkage
vr_shrinkage <- function(matU,
                         posU = matU > 0,
                         exclude_row = NULL,
                         exclude_col = NULL,
                         weights_col = NULL,
                         surv_only_na = TRUE) {
  
  vr_vec <- vr_vec_shrinkage(matU = matU,
                             posU = posU,
                             exclude_row = exclude_row,
                             exclude_col = exclude_col,
                             surv_only_na = surv_only_na)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr
#' @export vr_stasis
vr_stasis <- function(matU,
                      posU = matU > 0,
                      exclude_col = NULL,
                      weights_col = NULL,
                      surv_only_na = TRUE) {
  
  vr_vec <- vr_vec_stasis(matU = matU,
                          posU = posU,
                          exclude_col = exclude_col,
                          surv_only_na = surv_only_na)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr
#' @export vr_dorm_enter
vr_dorm_enter <- function(matU,
                          posU = matU > 0,
                          dorm_stages,
                          weights_col = NULL) {
  
  vr_vec <- vr_vec_dorm_enter(matU = matU,
                              posU = posU,
                              dorm_stages = dorm_stages)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr
#' @export vr_dorm_exit
vr_dorm_exit <- function(matU,
                         posU = matU > 0,
                         dorm_stages,
                         weights_col = NULL) {
  
  vr_vec <- vr_vec_dorm_exit(matU = matU,
                             posU = posU,
                             dorm_stages = dorm_stages)
  
  return(column_weight(vr_vec, weights_col))
}


#' @rdname vr
#' @export vr_fecundity
vr_fecundity <- function(matU,
                         matR,
                         posR = matR > 0,
                         exclude_col = NULL,
                         weights_row = NULL,
                         weights_col = NULL) {
  
  vr_vec <- vr_vec_fecundity(matU = matU,
                             matR = matR,
                             posR = posR,
                             exclude_col = exclude_col,
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
  if (is.nan(out)) out <- NA_real_
  return(out)
}
