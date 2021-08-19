#' @title Convert matrix population model into U, F and C matrices
#'
#' @description
#' Splits a matrix population model into three constituent matrices, \bold{U}
#' (growth and survival processes), \bold{F} (sexual reproduction) and \bold{C}
#' (clonal reproduction). \strong{Warning!} The functionality is very basic: it
#' assumes that sexual reproduction is located in the top row of the matrix, and
#' that everything else is growth or survival (i.e. the \bold{U} matrix).
#' Clonality is assumed to be non-existent.
#'
#' @param matA A matrix population model (i.e., a square projection matrix).
#' 
#' @return A list of three matrices: \code{matU},\code{matF} and \code{matC}. 
#' \code{matC} will always contain only zeros.
#' 
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' 
#' @family transformation
#' 
#' @examples
#' matA <- rbind(c(0.1,   0, 5.3, 4.2),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#' 
#' mpm_split(matA)
#' 
#' @export mpm_split
mpm_split <- function(matA) {
  
  # validate arguments
  checkValidMat(matA, fail_all_na = FALSE, fail_any_na = FALSE,
                warn_all_zero = FALSE)
  
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
  
  # copy attributes of matA
  attributes(matU) <- attributes(matA)
  attributes(matF) <- attributes(matA)
  attributes(matC) <- attributes(matA)
  
  return(list(matU = matU, matF = matF, matC = matC))
}
