#' Perturbation analysis of vital rates in a matrix population model
#'
#' @description 
#' Perturbs lower-level vital rates within a matrix population model and
#' measures the response (sensitivity or elasticity) of the per-capita 
#' population growth rate at equilibrium (\eqn{\lambda}), or, with a 
#' user-supplied function, any other demographic statistic.
#' 
#' These decompositions assume that all transition rates are products of a
#' stage-specific survival term (column sums of \code{matU}) and a lower level
#' vital rate that is conditional on survival (growth, shrinkage, stasis,
#' dormancy, or reproduction). Reproductive vital rates that are not conditional
#' on survival (i.e., within a stage class from which there is no survival) are
#' also allowed.
#'
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g., 
#'   progression, stasis, and retrogression).
#' @param matF The sexual component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions due to sexual reproduction).
#' @param matC The clonal component of a matrix population model (i.e., a square
#'   projection matrix reflecting transitions due to clonal reproduction).
#'   Defaults to \code{NULL}, indicating no clonal reproduction (i.e., 
#'   \code{matC} is a matrix of zeros).
#' @param pert Magnitude of the perturbation. Defaults to \code{1e-6}.
#' @param type Whether to return \code{sensitivity} or \code{elasticity} values. Defaults
#'   to \code{sensitivity}.
#' @param demog_stat The demographic statistic to be used, as in "the
#'   sensitivity/elasticity of \code{demog_stat} to vital rate perturbations."
#'   Defaults to the per-capita population growth rate at equilibrium
#'   (\eqn{\lambda}). Also accepts a user-supplied function that performs a
#'   calculation on a projection matrix and returns a single numeric value.
#' @param ... Additional arguments passed to the function \code{demog_stat}.
#' 
#' @return A list with 5 elements:
#' \item{survival}{sensitivity or elasticity of \code{demog_stat} to survival}
#' \item{growth}{sensitivity or elasticity of \code{demog_stat} to growth}
#' \item{shrinkage}{sensitivity or elasticity of \code{demog_stat} to shrinkage}
#' \item{fecundity}{sensitivity or elasticity of \code{demog_stat} to sexual
#' fecundity}
#' \item{clonality}{sensitivity or elasticity of \code{demog_stat} to clonality}
#' 
#' @author Rob Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @family perturbation analysis
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
#' perturb_vr(matU, matF)
#' 
#' # use elasticities rather than sensitivities
#' perturb_vr(matU, matF, type = "elasticity")
#' 
#' # use a larger perturbation than the default
#' perturb_vr(matU, matF, pert = 0.01)
#' 
#' # calculate the sensitivity/elasticity of the damping ratio to vital rate
#' #  perturbations
#' damping <- function(matA) {  # define function for damping ratio
#'   eig <- eigen(matA)$values
#'   dm <- rle(Mod(eig))$values
#'   return(dm[1] / dm[2])
#' }
#' 
#' perturb_vr(matU, matF, demog_stat = "damping")
#' 
#' @export perturb_vr
perturb_vr <- function(matU, matF, matC = NULL,
                       pert = 1e-6, type = "sensitivity",
                       demog_stat = "lambda", ...) {
  
  # validate arguments
  checkValidMat(matU)
  checkValidMat(matF)
  if (!is.null(matC)) checkValidMat(matC, warn_all_zero = FALSE)
  type <- match.arg(type, c("sensitivity", "elasticity"))
  
  # get statfun
  if (is.character(demog_stat) && demog_stat == "lambda") {
    statfun <- lambda
  } else {
    statfun <- try(match.fun(demog_stat), silent = TRUE)
    if (class(statfun) == "try-error") {
      stop(strwrap(prefix = " ", initial = "", "`demog_stat` must be `lambda` or the name 
                   of a function that returns a single numeric value.\n"), call. = FALSE)
    }
  }
  
  # matrix dimension
  m <- nrow(matU)
  
  # if matC null, convert to zeros
  if (is.null(matC)) matC <- matrix(0, nrow = m, ncol = m)
  
  # combine components into matA 
  matA <- matU + matF + matC
  
  # get sensitivity matrix
  sensA <- perturb_matrix(matA = matA, pert = pert, type = "sensitivity",
                          demog_stat = statfun, ...)
  
  stat <- statfun(matA, ...)
  
  # stage-specific survival
  sigma <- colSums(matU)
  
  # which stages defined only by age (single non-zero element in column of matU)
  stage_agedef <- apply(matU, 2, function(x) length(which(x > 0)) == 1)
  
  # survival-independent matA (divide columns of matA by corresponding element
  #  of sigma; retain original column of matA if no survival)
  noSurvA <- noSurvA_F <- t(t(matA) / sigma)
  noSurvA[,sigma == 0] <- matA[,sigma == 0]
  
  # sensitivity of survival
  matSensSurv <- noSurvA * sensA
  sensSurv <- colSums(matSensSurv)
  sensSurv[sigma == 0] <- 0  # zero out stages with no survival
  
  # lower and upper triangles (reflecting growth and retrogression)
  lwr <- upr <- matrix(0, nrow = m, ncol = m)
  lwr[lower.tri(lwr, diag = TRUE)] <- 1
  upr[upper.tri(upr, diag = TRUE)] <- 1
  
  # sensitivity of survival-independent U components
  matSensGrowShri <- matrix(NA_real_, m, m)
  
  for (i in 1:m) { 
    matSensGrowShri[,i] <- sensA[i,i] * (-sigma[i]) +
      sensA[,i] * (sigma[i])
  }
  
  # zero out non-existent U transitions, and age-defined stages
  matSensGrowShri[which(matU == 0)] <- 0
  matSensGrowShri[,stage_agedef == TRUE] <- 0
  
  # sensitivity to growth
  matSensGrow <- lwr * matSensGrowShri
  
  # sensitivity to shrinkage
  matSensShri <- upr * matSensGrowShri
  
  # modify sigma (column-specific survival) for reproduction transitions;
  # convert 0s to 1s, so as not to 'pull out' survival from columns
  #  from which there was no survival
  sigma_rep <- sigma
  sigma_rep[sigma_rep == 0] <- 1
  
  # sensitivity to fecundity
  matSensFec <- sensA
  matSensFec[which(matF == 0)] <- 0 # zero out non-existent transitions
  matSensFec <- t(t(matSensFec) * sigma_rep) # multiply columns by sigma_rep 
  
  # sensitivity to clonality
  matSensClo <- sensA
  matSensClo[which(matC == 0)] <- 0 # zero out non-existent transitions
  matSensClo <- t(t(matSensClo) * sigma_rep) # multiply columns by sigma_rep 
  
  if (type == "sensitivity") {
    
    out <- list(survival = sum(sensSurv, na.rm = TRUE),
                growth = sum(matSensGrow, na.rm = TRUE),
                shrinkage = sum(matSensShri, na.rm = TRUE),
                fecundity = sum(matSensFec, na.rm = TRUE),
                clonality = sum(matSensClo, na.rm = TRUE))
    
  } else {
    
    # elasticity of survival
    elasSurv <- sigma * sensSurv / stat
    
    # survival-independent U elasticities
    matElasGrowShri <- noSurvA * matSensGrowShri / stat
    
    # elasticity to growth and shrinkage
    matElasGrow <- lwr * matElasGrowShri
    matElasShri <- upr * matElasGrowShri
    
    # zero out grow/shri elasticities for stages with no survival
    # PB: I think this is unnecessary... just being 'safe'
    # elasGrow[sigma == 0] <- 0
    # elasShri[sigma == 0] <- 0
    
    # elasticity to fecundity and clonality
    matElasFec <- noSurvA * matSensFec / stat
    matElasClo <- noSurvA * matSensClo / stat
    
    out <- list(survival = sum(elasSurv, na.rm = TRUE),
                growth = sum(matElasGrow, na.rm = TRUE),
                shrinkage = sum(matElasShri, na.rm = TRUE),
                fecundity = sum(matElasFec, na.rm = TRUE),
                clonality = sum(matElasClo, na.rm = TRUE))
  }
  
  return(out)
}
