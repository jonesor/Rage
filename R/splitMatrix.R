#' Splits a matrix population model into the constituent U, F and C matrices
#'
#' This function splits a matrix population model into three constituent
#' matrices, U (growth and survival processes), F (sexual reproduction) and C
#' (clonal reproduction). Warning! The functionality is very basic — it assumes
#' that sexual reproduction is located in the top row of the matrix, and that
#' everything else is growth or survival (the U matrix). Clonality is assumed to
#' be non-existant.
#'
#' @param matA A matrix population model
#' @return A list of three matrices: \code{matU},\code{matF} and \code{matC}
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @examples
#' matA <- rbind(c(0.0, 0.0, 0.2, 1.3),
#'               c(0.5, 0.0, 0.0, 0.0),
#'               c(0.0, 0.3, 0.0, 0.0),
#'               c(0.0, 0.0, 0.1, 0.1))
#' 
#' splitMatrix(matA)
#' @export splitMatrix
splitMatrix <- function(matA) {
  
  # matrix dimension
  dim <- nrow(matA)
  
  # matU takes everything from matA except first row
  matU <- matA
  matU[1,] <- rep(0, dim)
  
  # matF takes first row of matA
  matF <- matrix(0, nrow = dim, ncol = dim)
  matF[1,] <- matA[1,]
  
  # matC is all zeros
  matC <- matrix(0, nrow = dim, ncol = dim)
  
  return(list(matU, matF, matC))
}
