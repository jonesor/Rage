#' Perturbation analysis of vital rates in a matrix population model
#'
#' Perturbs lower-level vital rates within a matrix population model and
#' measures the response of the per-capita population growth rate at equilibrium
#' (\eqn{\lambda}), or, with a user-supplied function, any other demographic
#' statistic.
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matF The sexual component of a matrix population model (i.e. a square
#'   projection matrix reflecting transitions due to sexual reproduction)
#' @param matC The clonal component of a matrix population model (i.e. a square
#'   projection matrix reflecting transitions due to clonal reproduction).
#'   Defaults to \code{NULL}, indicating no clonal reproduction (i.e.
#'   \code{matC} is a matrix of zeros).
#' @param pert Magnitude of the perturbation (defaults to \code{1e-6}).
#' @param demogstat The demographic statistic to be used, as in "the
#'   sensitivity/elasticity of ___ to vital rate perturbations." Defaults to the
#'   per-capita population growth rate at equilibrium (\eqn{lambda}). Also
#'   accepts a user-supplied function that performs a calculation on a
#'   projection matrix and returns a single numeric value.
#' @param ... Additional arguments passed to the function \code{demogstat}.
#' @return A data frame with 1 row and 10 columns:
#' \item{SSurvival}{sensitivity of \code{demogstat} to survival}
#' \item{SGrowth}{sensitivity of \code{demogstat} to growth}
#' \item{SShrinkage}{sensitivity of \code{demogstat} to shrinkage}
#' \item{SReproduction}{sensitivity of \code{demogstat} to sexual reproduction}
#' \item{SClonality}{sensitivity of \code{demogstat} to clonality}
#' \item{ESurvival}{elasticity of \code{demogstat} to survival}
#' \item{EGrowth}{elasticity of \code{demogstat} to growth}
#' \item{EShrinkage}{elasticity of \code{demogstat} to shrinkage}
#' \item{EReproduction}{elasticity of \code{demogstat} to sexual reproduction}
#' \item{EClonality}{elasticity of \code{demogstat} to clonality}
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
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
#' vitalRatePerturbation(matU, matF)
#' 
#' # use a larger perturbation than the default
#' vitalRatePerturbation(matU, matF, pert = 0.01)
#' 
#' # calculate the sensitivity/elasticity of the damping ratio to vital rate
#' #  perturbations
#' damping <- function(matA) {  # define function for damping ratio
#'   eig <- eigen(matA)$values
#'   dm <- rle(Mod(eig))$values
#'   return(dm[1] / dm[2])
#' }
#' 
#' vitalRatePerturbation(matU, matF, demogstat = 'damping')
#' 
#' @importFrom popbio lambda
#' @export vitalRatePerturbation
vitalRatePerturbation <- function(matU, matF, matC = NULL, pert = 1e-6,
                                  demogstat = "lambda", ...) {
  
  # get statfun
  if (is.function(demogstat)) {
    statfun <- demogstat
  } else if (demogstat == "lambda") {
    statfun <- popbio::lambda
  } else if (is.character(demogstat)) {
    statfun <- match.fun(demogstat, descend = FALSE)
  } else {
    stop("demogstat must be 'lambda' or the name of a function that returns a
          single numeric value")
  }
  
  # matrix dimension
  matDim <- nrow(matU)
  
  # if matC null, convert to zeros
  if (is.null(matC)) {
    matC <- matrix(0, matDim, matDim)
  }
  
  # combine components into matA 
  matA <- matU + matF + matC
  
  # statfun
  stat <- statfun(matA, ...)
  
  # initialize sensitivity matrix
  sensA <- matrix(NA, matDim, matDim)
  
  # matrix perturbation
  for (i in 1:matDim) {
    for (j in 1:matDim) {
      pertA <- matA
      pertA[i, j] <- pertA[i, j] + pert
      statPert <- statfun(pertA, ...)
      sensA[i, j] <- (stat - statPert) / (matA[i, j] - pertA[i, j])
    }
  }
  
  # stage-specific survival
  sigma <- colSums(matU)
  
  # which stages defined only by age (single non-zero element in column of matU)
  stage_agedef <- apply(matU, 2, function(x) length(which(x > 0))) == 1
  
  # survival-independent matA (divide columns of matA by corresponding element
  #  of sigma; retain original column of matA if no survival)
  noSurvA <- noSurvA_F <- t(t(matA) / sigma)
  noSurvA[,sigma == 0] <- matA[,sigma == 0]
  
  # sensitivity and elasticity of survival
  matSensSurv <- noSurvA * sensA
  sensSurv <- colSums(matSensSurv)
  sensSurv[sigma == 0] <- 0  # zero out stages with no survival
  elasSurv <- sigma * sensSurv / stat
  
  # lower and upper triangles (reflecting growth and retrogression)
  lwr <- upr <- matrix(0, nrow = matDim, ncol = matDim)
  lwr[lower.tri(lwr, diag = TRUE)] <- 1
  upr[upper.tri(upr, diag = TRUE)] <- 1
  
  # sensitivity of survival-independent U components
  matSensGrowShri <- matrix(NA_real_, matDim, matDim)
  
  for (i in 1:matDim) { 
    matSensGrowShri[,i] <- sensA[i,i] * (-sigma[i]) +
      sensA[,i] * (sigma[i])
  }
  
  # zero out non-existent U transitions, and age-defined stages
  matSensGrowShri[which(matU == 0)] <- 0
  matSensGrowShri[,stage_agedef == TRUE] <- 0
  
  # sensitivity to growth
  matSensGrow <- lwr * matSensGrowShri
  sensGrow <- colSums(matSensGrow, na.rm = TRUE)
  
  # sensitivity to shrinkage
  matSensShri <- upr * matSensGrowShri
  sensShri <- colSums(matSensShri, na.rm = TRUE)
  
  # modify sigma (column-specific survival) for reproduction transitions;
  # convert 0s to 1s, because don't want to 'pull out' survival from columns
  #  from which there was no survival
  sigma_rep <- sigma
  sigma_rep[sigma_rep == 0] <- 1
  
  # sensitivity to fecundity
  matSensFec <- sensA
  matSensFec[which(matF == 0)] <- 0 # zero out non-existent transitions
  matSensFec <- t(t(matSensFec) * sigma_rep) # multiply columns by sigma_rep 
  sensFec <- colSums(matSensFec)
  
  # sensitivity to clonality
  matSensClo <- sensA
  matSensClo[which(matC == 0)] <- 0 # zero out non-existent transitions
  matSensClo <- t(t(matSensClo) * sigma_rep) # multiply columns by sigma_rep 
  sensClo <- colSums(matSensClo)
  
  # survival-independent U elasticities
  matElasGrowShri <- noSurvA * matSensGrowShri / stat
  
  # elasticity to growth
  matElasGrow <- lwr * matElasGrowShri
  elasGrow <- colSums(matElasGrow, na.rm = TRUE)
  
  # elasticity to shrinkage
  matElasShri <- upr * matElasGrowShri
  elasShri <- colSums(matElasShri, na.rm = TRUE)
  
  # zero out grow/shri elasticities for stages with no survival
  # PB: I think this is unnecessary... just being 'safe'
  elasGrow[sigma == 0] <- 0
  elasShri[sigma == 0] <- 0
  
  # elasticity to fecundity
  matElasFec <- noSurvA * matSensFec / stat
  elasFec <- colSums(matElasFec)
  
  # elasticity to clonality
  matElasClo <- noSurvA * matSensClo / stat
  elasClo <- colSums(matElasClo)
  
  out <- data.frame("SSurvival" = sum(sensSurv, na.rm = TRUE),
                    "SGrowth" = sum(sensGrow, na.rm = TRUE),
                    "SShrinkage" = sum(sensShri, na.rm = TRUE),
                    "SReproduction" = sum(sensFec, na.rm = TRUE),
                    "SClonality" = sum(sensClo, na.rm = TRUE),
                    "ESurvival" = sum(elasSurv, na.rm = TRUE),
                    "EGrowth" = sum(elasGrow, na.rm = TRUE),
                    "EShrinkage" = sum(elasShri, na.rm = TRUE),
                    "EReproduction" = sum(elasFec, na.rm = TRUE),
                    "EClonality" = sum(elasClo, na.rm = TRUE))
  
  return(out)
}
