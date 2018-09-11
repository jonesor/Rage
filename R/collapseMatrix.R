#' Collapse a matrix population model to a smaller number of stages
#'
#' Collapse a matrix population model to a smaller number of stages. For
#' instance, to compare properties of multiple projection matrices with
#' different numbers of stages, one might first collapse those matrices to a
#' standardized set of stages (e.g. propagule, pre-reproductive, reproductive,
#' and post-reproductive). The transition rates in the collapsed matrix are a
#' weighted average of the transition rates from the relevant stages of the
#' original matrix, weighted by the relative proportion of each stage class
#' expected at the stable distribution.
#' 
#' @param matU A square matrix containing only survival-related transitions
#'   (i.e. progression, stasis, retrogression).
#' @param matF A square matrix containing only sexual reproduction-related
#'   transitions.
#' @param matC A square matrix containing only clonal reproduction-related
#'   transitions. Defaults to \code{NULL}, indicating no clonal reproduction
#'   (i.e. \code{matC} is a matrix of zeros).
#' @param collapse A list giving the mapping between stages of the original
#'   matrix and the desired stages of the collapsed matrix (e.g. \code{list(1,
#'   2:3, 4)}). The indices of \code{collapse} correspond to the desired stages
#'   of the collapsed matrix, and the corresponding values give the stage index
#'   or vector of stage indices from the original matrix that correspond to the
#'   relevant stage of the collapsed matrix.
#' @return A list with four elements:
#'   \item{matA}{Collapsed projection matrix}
#'   \item{matU}{Survival component of the collapsed projection matrix}
#'   \item{matF}{Sexual reproduction component of the collapsed projection matrix}
#'   \item{matC}{Clonal reproduction component of the collapsed projection matrix}
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @references Salguero-Gomez, R. & Plotkin, J. B. (2010) Matrix dimensions bias
#'   demographic inferences: implications for comparative plant demography. The
#'   American Naturalist 176, 710-722.
#' @note This method of collapsing a matrix population model preserves the
#'   equilibrium population growth rate (\eqn{lamda}) and relative stable
#'   distribution, but is not expected to preserve other traits such as relative
#'   reproductive values, sensitivities, net reproductive rates, life
#'   expectancy, etc.
#' @seealso \code{\link{standardizeMatrix}}
#' @examples
#' matU <- rbind(c(  0,   0,    0,    0),
#'               c(0.5,   0,    0,    0),
#'               c(  0, 0.3,    0,    0),
#'               c(  0,   0,  0.2,  0.1))
#' 
#' matF <- rbind(c(  0,   0,  1.1,  1.6),
#'               c(  0,   0,  0.8,  0.4),
#'               c(  0,   0,    0,    0),
#'               c(  0,   0,    0,    0))
#'               
#' matC <- rbind(c(  0,   0,  0.4,  0.5),
#'               c(  0,   0,  0.3,  0.1),
#'               c(  0,   0,    0,    0),
#'               c(  0,   0,    0,    0))
#' 
#' # collapse reproductive stages
#' collapse1 <- list(1, 2, 3:4)
#' collapseMatrix(matU, matF, matC, collapse1)
#' 
#' # collapse pre-reproductive stages, and reproductive stages
#' collapse2 <- list(1:2, 3:4)
#' collapseMatrix(matU, matF, matC, collapse2)
#' @export collapseMatrix
collapseMatrix <- function(matU, matF, matC = NULL, collapse) {
  
  # populate matC with zeros, if NULL
  if (is.null(matC)) {
    matC <- matrix(0, nrow = nrow(matU), ncol = ncol(matU))
  }
  
  # sum components to matA
  matA <- matU + matF + matC
  
  # ensure no NA
  if (any(is.na(matA))) {
    stop("Cannot collapse projection matrix containing NAs")
  }
  
  # dimensions of original and collapse matrices
  originalDim <- nrow(matA)
  collapseDim <- length(collapse)
  
  P <- matrix(0, nrow = collapseDim , ncol = originalDim)
  
  for (i in 1:collapseDim) {
    columns <- as.numeric(collapse[[i]])
    if (!is.na(columns[1])) {
      P[i, columns] <- 1
    }
  }
  
  Q <- t(P)
  w <- Re(eigen(matA)$vectors[ ,which.max(Re(eigen(matA)$values))])
  w <- w / sum(w)
  
  columns <- which(colSums(Q) > 1)
  
  for (j in columns) {
    rows <- which(Q[,j] == 1)
    for (i in rows) {
      Q[i,j] <- w[i] / sum(w[rows])
    }
  }
  
  # collapse
  collapseA <- P %*% matA %*% Q
  collapseU <- P %*% matU %*% Q
  collapseF <- P %*% matF %*% Q
  collapseC <- P %*% matC %*% Q
  
  return(list(matA = collapseA, 
              matU = collapseU, 
              matF = collapseF,
              matC = collapseC))
}
