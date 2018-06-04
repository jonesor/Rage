#' Calculates Keyfitz' entropy
#' 
#' This function calculates Keyfitz' entropy from an lx
#' (survivorship) vector with even intervals derived from a matrix population model.
#' 
#' @param matU A matrix containing only survival-dependent processes (e.g. progression,
#' stasis, retrogression).
#' @param startLife The first stage at which the author considers the beginning
#' of life in the life cycle of the species. It defaults to the first stage.
#' @param nSteps A cutoff for the decomposition of age-specific survival ('lx'), and when pertinent,
#' for age-specific sexual reproduction ('mx') and age-specific clonal reproduction ('cx'). This allows
#' excluding mortality and fertility plateaus. See function 'qsdConverge' for more information. When not
#' specified, this argument defaults to 100.
#' @param trapeze A logical argument indicating whether the trapezoidal
#' approximation should be used for approximating the definite integral.
#' @param trimlx A numerical argument determining how lx should be trimmed. 
#' For example, a value of 0.05 means that lx values are retained only if they 
#' are greater than or equal to 0.05. NB: 0 values in lx mean that Keyfitz 
#' entropy cannot be estimated.
#' @return Returns an estimate of Keyfitz' life table entropy based on an lx
#' (survivorship) vector obtained from matU
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' @references  %% ~~references~~
#' @examples
#'
#' matU <- matrix (c(0, 0, 0, 0, 0.6, 0, 0, 0, 0, 0.4, 0, 0, 0, 0, 0.7, 0.1), nrow = 4, byrow = TRUE)
#' kEntropy(matU, nSteps=100)
#' kEntropy(matU,trapeze=FALSE)
#' 
#' matU <- matrix (c(0.2, 0, 0, 0, 0.3, 0.4, 0.1, 0, 0.1, 0.1, 0.2, 0.3, 0, 0.2, 0.6, 0.5), nrow = 4, byrow = TRUE)
#' kEntropy(matU, nSteps = 10)
#' kEntropy(matU, nSteps = 20)
#' kEntropy(matU, nSteps = 100)
#' kEntropy(matU, nSteps = 100, trapeze=TRUE)
#' 
#' @export kEntropy
kEntropy <- function(matU, startLife = 1, nSteps = 1000, trapeze = FALSE, trimlx = 0.05){
  
  if (dim(matU)[1]!=dim(matU)[2]) stop("Your matrix population model is not a square matrix")
  if (any(is.na(matU))) stop("NAs exist in matU")
  if (length(which(colSums(matU)>1))>0) print("Warning: matU has at least one stage-specific survival value > 1")
  
  #Age-specific survivorship (lx) (See top function on page 120 in Caswell 2001):
  matDim = dim(matU)[1]
  matUtemp = matU
  survivorship = array(NA, dim = c(nSteps, matDim))
  for (o in 1:nSteps){
    survivorship[o, ] = colSums(matUtemp %*% matU)
    matUtemp = matUtemp %*% matU
  }
  
  lx = survivorship[, startLife]
  lx = c(1, lx[1:(length(lx) - 1)])
  
  #Trim lx
  lx <- lx[lx>=trimlx]
  
  if(trapeze == TRUE){
    ma <- function(x,n=2){filter(x,rep(1/n,n), sides=2)}
    lx2 <- na.omit(as.vector(ma(lx)))
    return(-sum(lx2*log(lx2))/sum(lx2))
    }else{
      return(-sum(lx*log(lx))/sum(lx))
    }
  }
