#' A function to perform perturbation of vital rates of the matrix model for any demographic statistic.
#' 
#' A function to perform perturbation of vital rates of the matrix model and
#' measure the response of the per-capita population growth rate at equilibrium
#' or (with a user-supplied function) any other demographic statistic.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param matU The U matrix (processes related to survival, growth and retrogression).
#' @param matF The F matrix (sexual reproduction processes).
#' @param matC The C matrix (clonal reproduction processes).
#' @param demogstat The demographic statistic to be used, as
#'  in "the sensitvity/elasticity of ___ to vital rate perturbations."
#'  Defaults to the per-capita population growth rate at equilibrium. 
#'  Also accepts a user-supplied function that performs a calculation on
#'  a projection matrix and returns a single numeric value.
#' @param pert Magnitude of the perturbation (defaults to 0.001).
#' @return A data frame containing...
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gomez <r.salguero@@sheffield.ac.uk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' \dontrun{
#' data(Compadre)
#' pira <- subsetDB(Compadre, SpeciesAccepted == "Pinus radiata")
#' 
#' vitalRatePerturbation(pira)  # Not implemented yet!
#' 
#' vitalRatePerturbation(matU = pira@matU, matF = pira@matF)
#' vitalRatePerturbation(matU = pira@matU, matF = pira@matF, pert = 1e-03)
#' # use a smaller perturbation
#' vitalRatePerturbation(matU = pira@matU, matF = pira@matF, pert = 1e-10)
#' 
#' # calculate the sensitivity/elasticity of the damping ratio to vital rate perturbations
#' damping <- function(matU, matF, matC){  # define function for damping ratio
#' A <- matU+matF+matC
#' eig <- eigen(A)$values
#' dm <- rle(Mod(eig))$values
#' return(dm[1] / dm[2])
#' }
#' vitalRatePerturbation(matU = pira@matU, matF = pira@matF, stat = damping)
#' vitalRatePerturbation(matU = pira@matU, matF = pira@matF, stat = "damping")
#' 
#' }
#' 
#' 
#' @importFrom popbio eigen.analysis
#' @export vitalRatePerturbation
vitalRatePerturbation <- function(matU, matF, matC = NULL, demogstat = "lambda",
                                  pert = 1e-03){
  #Function to calculate vital rate level sensitivities and elasticities
  
  # If matC is actually NULL, then then the matA calculation returns
  # integer(0), and the whole thing fails. Thus, coerce to matrix
  # of 0s of same dimension as matU. Kind of a hacky fix, but it 
  # seems to work.
  if(is.null(matC)){
    matC <- matrix(0, nrow(matU), dim(matA)[1])
  }
  matA <- matU + matF + matC
  aDim <- dim(matA)[1]
  fakeA <- matA
  sensA <- elasA <- matrix(NA, aDim, aDim)
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
  propU[is.nan(propU)] <- 0
  propProg <- propRetrog <- propU
  propProg[upper.tri(propU, diag = TRUE)] <- 0
  propRetrog[lower.tri(propU, diag = TRUE)] <- 0
  propStasis <- matrix(diag(aDim) * diag(propU), aDim, aDim)
  propF <- matF / matA
  propF[is.nan(propF)] <- 0
  propC <- matC / matA
  propC[is.nan(propC)] <- 0
  
  for (i in 1:aDim){
    for (j in 1:aDim){
      fakeA <- matA
      fakeA[i, j] <- fakeA[i, j] + pert
      statPert <- statfun(fakeA, ...)
      sensA[i, j] <- (stat - statPert) / (matA[i, j] - fakeA[i, j])
    }
  }
  
  sensA <- Re(sensA)
  
  #Survival-independent A matrix
  uIndep <- matrix(NA, aDim, aDim)
  u <- colSums(matU)
  for (j in which(u > 0)) uIndep[, j] <- matA[, j] / u[j]
  
  sensSigmaA <- uIndep * sensA
  
  #Little fix for semelparous species
  uPrime <- u
  #uPrime[u==0] <- 0.001
  elasSigmaA <- t(t(sensSigmaA) * uPrime) / stat
  
  elasA <- sensA * matA / stat
  
  #Extracting survival vital rate
  uDistrib <- matrix(0, ncol = aDim, nrow = aDim)
  for (j in which(u > 0)) uDistrib[, j] <- matU[, j] / u[j]
  #Extracting fecundity vital rates:
  f <- colSums(matF)
  fDistrib <- matrix(0, ncol = aDim, nrow = aDim)
  for (j in which(f > 0)) fDistrib[, j]=matF[, j] / f[j]
  #Extracting clonality vital rates:
  c <- colSums(matC)
  cDistrib <- matrix(0, ncol = aDim, nrow = aDim)
  for (j in which(c > 0)) cDistrib[, j] = matC[, j] / c[j]
  
  
  SuDistrib=sensA * uDistrib
  SfDistrib=sensA * fDistrib
  ScDistrib=sensA * cDistrib
  
  
  out <- data.frame(SSurvival = NA, SGrowth = NA, SShrinkage = NA, SReproduction = NA,
                    SClonality = NA, ESurvival = NA, EGrowth = NA, EShrinkage = NA,
                    EReproduction = NA, EClonality = NA)
  
  
  #Still to be done
  out$SSurvival <- sum(sensSigmaA, na.rm = TRUE)
  out$SGrowth <- sum(sensA * uDistrib * propProg, na.rm = TRUE)
  out$SShrinkage <- sum(sensA * uDistrib * propRetrog, na.rm = TRUE)
  out$SReproduction <- sum(sensA * fDistrib * propF, na.rm = TRUE)
  out$SClonality <- sum(sensA * cDistrib * propC, na.rm = TRUE)
  
  EuDistrib <- sensA * uDistrib * matrix(u, nrow = aDim, ncol = aDim, byrow = TRUE) / stat
  EfDistrib <- sensA * fDistrib * matrix(f, nrow = aDim, ncol = aDim, byrow = TRUE) / stat
  EcDistrib <- sensA * cDistrib * matrix(c, nrow = aDim, ncol = aDim, byrow = TRUE) / stat
  
  out$ESurvival <- sum(elasSigmaA, na.rm = TRUE)
  out$EGrowth <- sum(EuDistrib * propProg, na.rm = TRUE)
  out$EShrinkage <- sum(EuDistrib * propRetrog, na.rm = TRUE)
  out$EReproduction <- sum(EfDistrib * propF, na.rm = TRUE)
  out$EClonality <- sum(EcDistrib * propC, na.rm = TRUE)
  
  return(out) 
}
