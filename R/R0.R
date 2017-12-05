#' Calculates net reproductive value
#' 
#' Calculates net reproductive value from a matU
#' (survival-dependent processes) and either a matF (sexual reproduction) and/or a
#' matC (clonal reproduction).
#' 
#' @param matU A matrix containing only survival-dependent processes (e.g. progression,
#' stasis, retrogression).
#' @param matF A matrix containing only sexual reproduction, with zeros
#' elsewhere.
#' @param matC A matrix containing only clonal reproduction, with zeros
#' elsewhere. If not provided, it defaults to a matrix with all zeros.
#' @param startLife The first stage at which the author considers the beginning
#' of life in the life cycle of the species. It defaults to the first stage.
#' @return Returns the net reproductive value of the matrix population model. When both 'matF'
#' and 'matC' are provided, it outputs the net reprodutive value for sexual
#' reproduction only, for clonal reproduction only, and for both types of
#' reproduction together.
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' Hal Caswell <h.caswell@@uva.nl>
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#' Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#' 978-0878930968
#' @examples
#' 
#' matU <- matrix (c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
#' matF <- matrix (c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
#' matC <- matrix (c(0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
#' 
#' R0(matU, matF)
#' R0(matU, matC)
#' R0(matU, matF, matC, startLife=1)
#' R0(matU, matF, matU, startLife=4)
#' 
#' @import MASS
#' @export R0

R0 <- function(matU, matF, matC=FALSE, startLife=1){
  
  #Error checks
  if (dim(matU)[1]!=dim(matU)[2]) stop("Your matrix population model is not a square matrix")
  if (any(is.na(matF)) & any(is.na(matC))) stop("NAs exist in both matF and matC")
  if (sum(matC)==0) {matC=matrix(0,dim(matU),dim(matU))}
  if (length(which(colSums(matU)>1))>0) print("Warning: matU has at least one stage-specific survival value > 1")
  
  #The full matrix population model
  matA = matU + matF +matC
  
  R0=NULL
  matDim <- dim(matA)[1]
  N <- solve(diag(matDim)-matU)  #Fundamental matrix, which states the amount of time units spent in each stage on average
  
  if (sum(matF)>0){
    R0_matF <- matF%*%N
    R0$Fec <- R0_matF[startLife,startLife]
  }
  if (sum(matC)>0){
    R0_matC <- matC%*%N
    R0$Clo <- R0_matC[startLife,startLife]
  }
  if (sum(matF)>0 & sum(matC)>0){
    matFC <- matF + matC
    R0_matFC <- matFC%*%N
    R0$FecClo <- R0_matFC[startLife,startLife]
  }
  
return(R0)
}
