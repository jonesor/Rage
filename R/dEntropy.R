#' Calculates Demetrius' entropy
#' 
#' This function calculates Demetrius' entropy from a matrix population model
#' 
#' @param matU A matrix containing only survival-dependent processes (e.g. progression,
#' stasis, retrogression).
#' @param matF A matrix containing only sexual reproduction, with zeros
#' elsewhere.
#' @param matC A matrix containing only clonal reproduction, with zeros
#' elsewhere. If not provided, it defaults to a matrix with all zeros.
#' @param startLife The first stage at which the author considers the beginning
#' of life in the life cycle of the species. It defaults to the first stage.
#' @param nSteps A cutoff for the decomposition of age-specific survival ('lx'), and when pertinent,
#' for age-specific sexual reproduction ('mx') and age-specific clonal reproduction ('cx'). This allows
#' excluding mortality and fertility plateaus. See function 'qsdConverge' for more information. When not
#' specified, this argument defaults to 100.
#' @return Returns an estimate of Demetrius' entropy. When both 'matF' and 'matC' are provided,
#' it outputs Demetrius' entropy for sexual reproduction only ('Fec'), for clonal reproduction only ('Clo'),
#' and for both types of reproduction together ('FecClo').
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references L. Demetrius. 1978. Adaptive value, entropy and survivorship
#' curves. Nature 275, 213 - 214. doi:10.1038/275213a0
#' @examples
#' matU <- matrix (c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = TRUE)
#' matF <- matrix (c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = TRUE)
#' matC <- matrix (c(0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = TRUE)
#' 
#' dEntropy(matU, matF, matC, nSteps=5)
#' dEntropy(matU, matF, matC, nSteps=10)
#' dEntropy(matU, matF, matC, nSteps=100)
#' 
#' #' @import MASS
#' @export dEntropy
dEntropy <- function(matU, matF, matC=FALSE, startLife=1, nSteps=1000){
  
  #Error checks
  if (dim(matU)[1]!=dim(matU)[2]) stop("Your matrix population model is not a square matrix")
  if (any(is.na(matF)) & any(is.na(matC))) stop("NAs exist in both matF and matC")
  if (sum(matC)==0) {matC=matrix(0,dim(matU),dim(matU))}
  if (length(which(colSums(matU)>1))>0) print("Warning: matU has at least one stage-specific survival value > 1")
  
  matA = matU + matF + matC
  matDim = dim(matA)[1]
  
  dEntropy=NULL
  
  #population growth rate in log scale (rmax)
  r = log(max(Re(eigen(matA)$value)))  
  
  #Age-specific survivorship (lx) (See top function on page 120 in Caswell 2001):
  matUtemp = matU
  survivorship = array(NA, dim = c(nSteps, matDim))
  for (o in 1:nSteps){
    survivorship[o, ] = colSums(matUtemp %*% matU)
    matUtemp = matUtemp %*% matU
  }
  
  lx = survivorship[, startLife]
  lx = c(1, lx[1:(length(lx) - 1)])
  
  if(!missing(matF)){
    if(sum(matF,na.rm=TRUE)==0){
      warning("matF contains only 0 values")
    }
    #Age-specific fertility (mx, Caswell 2001, p. 120)
    ageFertility = array(0, dim = c(nSteps, matDim))
    fertMatrix = array(0, dim = c(nSteps, matDim))
    matUtemp2 = matU
    e = matrix(rep(1, matDim))
    for (q in 1:nSteps) {
      fertMatrix = matF %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
      ageFertility[q, ] = colSums(fertMatrix)
      matUtemp2 = matUtemp2 %*% matU
    }  
    mx = ageFertility[, startLife]
    mx = c(0, mx[1:(length(mx) - 1)])
    mx[is.nan(mx)]=0
  }
  
  if(!missing(matC)){
    if(sum(matC,na.rm=TRUE)==0){
      warning("matC contains only 0 values")
    }
    #Age-specific clonality (cx)
    ageClonality = array(0, dim = c(nSteps, matDim))
    clonMatrix = array(0, dim = c(nSteps, matDim))
    matUtemp2 = matU
    e = matrix(rep(1, matDim))
    for (q in 1:nSteps) {
      clonMatrix = matC %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
      ageClonality[q, ] = colSums(clonMatrix)
      matUtemp2 = matUtemp2 %*% matU
    }  
    cx = ageClonality[, startLife]
    cx = c(0, cx[1:(length(cx) - 1)])
    cx[is.nan(cx)]=0
  }
  
  if(!missing(matF) & !missing(matC)){
    mxcx=mx+cx
  }
  
  
  if (sum(mx)>0){
    limiteFx <- min(length(mx[which(!is.na(mx))]), length(lx[which(!is.na(lx))]))   #Calculating the last time interval at which both mx and lx were able to be calculated
    lxmx <- lx[1:limiteFx]*mx[1:limiteFx]  #Step-by-step multiplication of lx and mx
    lxmx[which(lxmx==0)] <- 1 #Coertion of the first few values for which mx = 0 to take values lxmx = 1, so that it can be log-scaled, below
    loglxmx <- log(lxmx)
    loglxmx[which(lxmx==0)] <- NA
    
    dEntropy$Fec <- abs(sum(lxmx*loglxmx)/sum(lxmx))
  }
  
  if (sum(cx)>0){
    limiteCx <- min(length(cx[which(!is.na(cx))]), length(lx[which(!is.na(lx))]))   #Calculating the last time interval at which both cx and lx were able to be calculated
    lxcx <- lx[1:limiteCx]*cx[1:limiteCx]  #Step-by-step multiplication of lx and cx
    lxcx[which(lxcx==0)] <- 1 #Coertion of the first few values for which cx = 0 to take values lxcx = 1, so that it can be log-scaled, below
    loglxcx <- log(lxcx)
    loglxcx[which(lxcx==0)] <- NA
    
    dEntropy$Clo <- abs(sum(lxcx*loglxcx)/sum(lxcx))
  }
  
  if (sum(mx)>0 & sum(cx)>0){
    limiteMxCx <- min(length(mx[which(!is.na(mx))]), length(cx[which(!is.na(cx))]), length(lx[which(!is.na(lx))]))   #Calculating the last time interval at which both mx, cx and lx were able to be calculated
    mxcx=rowSums(cbind(mx,cx))
    lxmxcx <- lx[1:limiteMxCx]*mxcx[1:limiteMxCx]  #Step-by-step multiplication of lx and cx
    lxmxcx[which(lxmxcx==0)] <- 1 #Coertion of the first few values for which cx = 0 to take values lxcx = 1, so that it can be log-scaled, below
    loglxmxcx <- log(lxmxcx)
    loglxmxcx[which(lxmxcx==0)] <- NA
    
    dEntropy$FecClo <- abs(sum(lxmxcx*loglxmxcx)/sum(lxmxcx))
  }
  return(dEntropy)
}
