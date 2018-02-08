#' A function to perform element perturbation of a matrix population model for 
#' any demographic statistic.
#' 
#' A function to perform element perturbation of a matrix population model and
#' measure the response of the per-capita population growth rate at equilibrium (λ)
#' or (with a user-supplied function) any other demographic statistic.
#' 
#' @param matU The U matrix (processes related to survival, growth and retrogression).
#' @param matF The F matrix (sexual reproduction processes).
#' @param matC The C matrix (clonal reproduction processes).
#' @param pert Perturbation parameter.
#' @param demogstat A character string that is the name of a function.
#' The default is the the per-capita population growth rate at equilibrium (λ). 
#'  Also accepts a user-supplied function that performs a calculation on
#'  a projection matrix and returns a single numeric value.
#' @param ... Additional arguments passed to the function \code{demogstat}.
#' @return %% ~Describe the value returned 
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gomez <r.salguero@@sheffield.ac.uk>
#' @references %% ~~references~~
#' @examples
#' 
#' 
#' @export matrixElementPerturbation
matrixElementPerturbation <- function(matU, matF, matC=NULL, demogstat = "lambda",
                                      pert = 0.001, ...){
  if(is.null(matC)){
    matC <- matrix(0, nrow(matU), ncol(matU))
  }
  matA <- matU + matF + matC
  aDim <- nrow(matA)
  fakeA <- matA
  sensA <- elasA <- matrix(NA,aDim,aDim)
  if (is.function(demogstat)) {
    statfun <- demogstat
  } else if (is.character(demogstat)) {
    statfun <- getFunction(demogstat)
  } else if (demogstat == "lambda") {
    statfun <- function(matA) Re(eigen(matA)$values[1])
  } else {
    stop("demogstat must be 'lambda' or the name of a function that returns a *single* numeric value.")
  }
  stat <- statfun(matA, ...)
  
  propU <- matU / matA
  propU[is.nan(propU)] <- NA
  propProg <- propRetrog <- propU
  propProg[upper.tri(propU, diag = TRUE)] <- NA
  propRetrog[lower.tri(propU, diag = TRUE)] <- NA
  propStasis <- matrix(diag(aDim) * diag(propU), aDim, aDim)
  propF <- matF / matA
  propF[is.nan(propF)] <- NA
  propC <- matC / matA
  propC[is.nan(propC)] <- NA
  
  for (i in 1:aDim){
    for (j in 1:aDim){
      fakeA <- matA
      fakeA[i, j] <- fakeA[i, j] + pert
      statPert <- statfun(fakeA, ...)
      sensA[i, j] <- (stat - statPert) / (matA[i, j] - fakeA[i, j])
    }
  }
  
  sensA <- Re(sensA)
  elasA <- sensA * matA / stat
  
  out <- data.frame(SStasis = NA, SProgression = NA, SRetrogression = NA,
                    SFecundity = NA, SClonality = NA, EStasis = NA,
                    EProgression = NA, ERetrogression = NA,
                    EFecundity = NA, EClonality = NA)
  
  out$SStasis <- sum(sensA * propStasis, na.rm = TRUE)
  out$SRetrogression <- sum(sensA * propRetrog, na.rm = TRUE)
  out$SProgression <- sum(sensA * propProg, na.rm = TRUE)
  out$SFecundity <- sum(sensA * propF, na.rm = TRUE)
  out$SClonality <- sum(sensA * propC, na.rm =  TRUE)
  out$EStasis <- sum(elasA * propStasis, na.rm = TRUE)
  out$EProgression <- sum(elasA * propProg, na.rm = TRUE)
  out$ERetrogression <- sum(elasA * propRetrog, na.rm = TRUE)
  out$EFecundity <- sum(elasA * propF, na.rm = TRUE)
  out$EClonality <- sum(elasA * propC, na.rm = TRUE)
  
  return(out) 
}
