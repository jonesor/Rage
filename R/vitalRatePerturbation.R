#' Perturbation analysis of vital rates in a matrix population model
#'
#' Perturbs lower-level vital rates within a matrix population model and
#' measures the response of the per-capita population growth rate at equilibrium
#' (\eqn{lambda}), or, with a user-supplied function, any other demographic
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
#'   per-capita population growth rate at equilibrium (\eqn{\lambda}). Also
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
#' vitalRatePerturbation(matU, matF, demogstat = damping)
#' 
#' @importFrom popbio lambda
#' @importFrom methods getFunction
#' @export vitalRatePerturbation
vitalRatePerturbation <- function(matU, matF, matC = NULL, pert = 1e-6,
                                  demogstat = "lambda", ...) {
  
  # get statfun
  if (is.function(demogstat)) {
    statfun <- demogstat
  } else if (demogstat == "lambda") {
    statfun <- popbio::lambda
  } else if (is.character(demogstat)) {
    statfun <- getFunction(demogstat)
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
  
  propU <- matU / matA
  propU[is.nan(propU)] <- 0
  propProg <- propRetrog <- propU
  propProg[upper.tri(propU, diag = TRUE)] <- 0
  propRetrog[lower.tri(propU, diag = TRUE)] <- 0
  propStasis <- diag(matDim) * diag(propU)
  propF <- matF / matA
  propF[is.nan(propF)] <- 0
  propC <- matC / matA
  propC[is.nan(propC)] <- 0
  
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
  
  # survival-independent A matrix
  uIndep <- matrix(NA, matDim, matDim)
  u <- colSums(matU)
  for(j in which(u > 0)) uIndep[, j] <- matA[, j] / u[j]
  
  sensSigmaA <- uIndep * sensA

  
  #Little fix for semelparous species
  # uPrime <- u
  #uPrime[u==0] <- 0.001
  elasSigmaA <- t(t(sensSigmaA) * u) / stat
  elasA <- sensA * matA / stat
  
  
  #Extracting survival vital rate
  uDistrib <- matrix(0, ncol = matDim, nrow = matDim)
  for (j in which(u > 0)) uDistrib[, j] <- matU[, j] / u[j]
  
  #Extracting fecundity vital rates:
  f <- colSums(matF)
  fDistrib <- matrix(0, ncol = matDim, nrow = matDim)
  for (j in which(f > 0)) fDistrib[, j]=matF[, j] / f[j]
  
  #Extracting clonality vital rates:
  c <- colSums(matC)
  cDistrib <- matrix(0, ncol = matDim, nrow = matDim)
  for (j in which(c > 0)) cDistrib[, j] = matC[, j] / c[j]

  
  SuDistrib <- sensA * uDistrib
  SfDistrib <- sensA * fDistrib
  ScDistrib <- sensA * cDistrib
  
  EuDistrib <- sensA * uDistrib * matrix(u, nrow = matDim, ncol = matDim, byrow = TRUE) / stat
  EfDistrib <- sensA * fDistrib * matrix(f, nrow = matDim, ncol = matDim, byrow = TRUE) / stat
  EcDistrib <- sensA * cDistrib * matrix(c, nrow = matDim, ncol = matDim, byrow = TRUE) / stat
  
  
  #Still to be done
  out <- data.frame(
    SSurvival = sum(sensSigmaA, na.rm = TRUE),
    SGrowth = sum(sensA * uDistrib * propProg, na.rm = TRUE),
    SShrinkage = sum(sensA * uDistrib * propRetrog, na.rm = TRUE),
    SReproduction = sum(sensA * fDistrib * propF, na.rm = TRUE),
    SClonality = sum(sensA * cDistrib * propC, na.rm = TRUE),
    ESurvival = sum(elasSigmaA, na.rm = TRUE),
    EGrowth = sum(EuDistrib * propProg, na.rm = TRUE),
    EShrinkage = sum(EuDistrib * propRetrog, na.rm = TRUE),
    EReproduction = sum(EfDistrib * propF, na.rm = TRUE),
    EClonality = sum(EcDistrib * propC, na.rm = TRUE)
  )
  
  return(out) 
}
