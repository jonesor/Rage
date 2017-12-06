#' Standardize vital rates
#'
#' @export
#' @param matU Survival matrix
#' @param matF Fecundity matrix
#' @param reproStages Logical vector indicating which stages reproductive
#' @param matrixStages Character vector of matrix stage types (e.g. "propagule",
#'   "active", or "dormant")
#' @return A list (TODO:)
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @examples
#' matU <- rbind( c(0, 0, 0, 0, 0), c(0.18, 0.16, 0, 0, 0), c(0.29, 0.23, 0.12,
#' 0, 0), c(0, 0, 0.34, 0.53, 0), c(0, 0, 0, 0.87, 0) )
#'
#' matF <- rbind( c(0, 0.13, 0, 0.96, 0), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0),
#' c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0) )
#'
#' reproStages <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
#' matrixStages <- c('prop', 'active', 'active', 'active', 'active')
#'
#' standardizedVitalRates(matU, matF, reproStages, matrixStages)
standardizedVitalRates <- function(matU, matF, reproStages, matrixStages) {

  # put non-reproductive stages at the end of the matrix
  rearr <- Rcompadre::rearrangeMatrix(matU, matF, reproStages, matrixStages)
  
  # defines which columns need to be collapsed for each of the four stages
  collapse <- reprodStages(rearr$matF,
                           rearr$nonRepInterRep,
                           rearr$reproStages,
                           rearr$matrixStages)
  
  # collapse
  xx <- Rcompadre::collapseMatrix(matU, matF, collapse = collapse)
  matUcollapse <- xx$matU
  matUcollapse[is.na(collapse), is.na(collapse)] <- NA
  
  matFcollapse <- xx$matF
  matFcollapse[is.na(collapse), is.na(collapse)] <- NA
  
  out <- vitalRates(matU = matUcollapse, matF = matFcollapse,
                    splitStages = rearr$matrixStages)
  return(out)
}
