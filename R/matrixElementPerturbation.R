#' Perturbation analysis of matrix elements in a matrix population model
#'
#' Perturbs elements within a matrix population model and measures the response
#' of the per-capita population growth rate at equilibrium (\eqn{\lambda}), or,
#' with a user-supplied function, any other demographic statistic.
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
#'   sensitivity/elasticity of ___ to matrix element perturbations." Defaults to
#'   the per-capita population growth rate at equilibrium (\eqn{lambda}). Also
#'   accepts a user-supplied function that performs a calculation on a
#'   projection matrix and returns a single numeric value.
#' @param ... Additional arguments passed to the function \code{demogstat}.
#' @return A data frame with 1 row and 10 columns:
#' \item{SStasis}{sensitivity of \code{demogstat} to stasis}
#' \item{SRetrogression}{sensitivity of \code{demogstat} to retrogression}
#' \item{SProgression}{sensitivity of \code{demogstat} to progression}
#' \item{SFecundity}{sensitivity of \code{demogstat} to sexual fecundity}
#' \item{SClonality}{sensitivity of \code{demogstat} to clonality}
#' \item{EStasis}{elasticity of \code{demogstat} to stasis}
#' \item{ERetrogression}{elasticity of \code{demogstat} to retrogression}
#' \item{EProgression}{elasticity of \code{demogstat} to progression}
#' \item{EFecundity}{elasticity of \code{demogstat} to sexual fecundity}
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
#' matrixElementPerturbation(matU, matF)
#' 
#' # use a larger perturbation than the default
#' matrixElementPerturbation(matU, matF, pert = 0.01)
#' 
#' # calculate the sensitivity/elasticity of the damping ratio to perturbations
#' damping <- function(matA) {  # define function for damping ratio
#'   eig <- eigen(matA)$values
#'   dm <- rle(Mod(eig))$values
#'   return(dm[1] / dm[2])
#' }
#' 
#' matrixElementPerturbation(matU, matF, demogstat = 'damping')
#' 
#' @importFrom popbio lambda
#' @export matrixElementPerturbation
matrixElementPerturbation <- function(matU, matF, matC = NULL, pert = 1e-6,
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
  
  propU <- matU / matA
  propU[is.nan(propU)] <- NA
  propProg <- propRetrog <- propU
  propProg[upper.tri(propU, diag = TRUE)] <- NA
  propRetrog[lower.tri(propU, diag = TRUE)] <- NA
  propStasis <- diag(matDim) * diag(propU)
  propF <- matF / matA
  propF[is.nan(propF)] <- NA
  propC <- matC / matA
  propC[is.nan(propC)] <- NA
  
  # initialize sensitivity matrix
  sensA <- matrix(NA, matDim, matDim)
  
  # matrix perturbation
  for (i in 1:matDim) {
    for (j in 1:matDim) {
      fakeA <- matA
      fakeA[i, j] <- fakeA[i, j] + pert
      statPert <- statfun(fakeA, ...)
      sensA[i, j] <- (stat - statPert) / (matA[i, j] - fakeA[i, j])
    }
  }
  
  sensA <- Re(sensA)
  elasA <- sensA * matA / stat
  
  out <- data.frame(
    SStasis = sum(sensA * propStasis, na.rm = TRUE),
    SRetrogression = sum(sensA * propRetrog, na.rm = TRUE),
    SProgression = sum(sensA * propProg, na.rm = TRUE),
    SFecundity = sum(sensA * propF, na.rm = TRUE),
    SClonality = sum(sensA * propC, na.rm =  TRUE),
    EStasis = sum(elasA * propStasis, na.rm = TRUE),
    EProgression = sum(elasA * propProg, na.rm = TRUE),
    ERetrogression = sum(elasA * propRetrog, na.rm = TRUE),
    EFecundity = sum(elasA * propF, na.rm = TRUE),
    EClonality = sum(elasA * propC, na.rm = TRUE)
  )
  
  return(out) 
}
