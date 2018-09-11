#' Generate a life table from a matrix population model
#'
#' This function uses age-from-stage decomposition methods to generate a life
#' table from a matrix population model. A detailed description of these methods
#' can be found in section 5.3 of Caswell (2001).
#'
#' @param matU A square matrix containing only survival-related transitions
#'   (i.e. progression, stasis, retrogression).
#' @param matF (Optional) A square matrix containing only sexual
#'   reproduction-related transitions.
#' @param matC (Optional) A square matrix containing only clonal
#'   reproduction-related transitions.
#' @param startLife The index of the first stage at which the author considers
#'   the beginning of life. Defaults to 1.
#' @param nSteps Number of time steps for which the life table will be
#'   calculated. Time steps are in the same units as the matrix population model
#'   (see MatrixPeriodicity column in metadata of COM(P)ADRE). Defaults to 1000.
#' @return A \code{data.frame} containing 7-13 columns and \code{nSteps} rows.
#'   Columns include:
#'   \item{x}{age}
#'   \item{lx}{survivorship to start of age class x}
#'   \item{dx}{proportion of original cohort dying in interval [x, x+1)}
#'   \item{qx}{force of mortality at age x}
#'   \item{Lx}{survivorship to middle of age class x}
#'   \item{Tx}{proportion of original cohort alive at age x and beyond}
#'   \item{ex}{remaining life expectancy at age x}
#'
#' If \code{matF} is provided, also includes:
#'   \item{mx}{per-capita rate of sexual reproduction at age x}
#'   \item{lxmx}{expected number of sexual offspring per original
#'   cohort member produced at age x}
#'   
#' If \code{matC} is provided, also includes:
#'   \item{cx}{per-capita rate of clonal reproduction at age x}
#'   \item{lxcx}{expected number of clonal offspring per original
#'   cohort member produced at age x}
#'   
#' If both \code{matF} and \code{matC} are provided, also includes:
#'   \item{mxcx}{per-capita rate of total reproduction (sexual + clonal) at age x}
#'   \item{lxmxcx}{expected number of total offspring (sexual + clonal) per original
#'   cohort member produced at age x}
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk> 
#' @author Hal Caswell <h.caswell@@uva.nl> 
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @references
#' Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and
#' Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968
#'
#' Caswell, H. (2006) Applications of Markov chains in demography. pp. 319-334
#' in A.N. Langville and W.J. Stewart (editors) MAM2006: Markov Anniversary
#' Meeting. Boson Books, Raleigh, North Caroline, USA
#'
#' Jones, O.R. et al. (2014) Diversity of ageing across the tree of life.
#' Nature, 505(7482), 169–173
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
#' matC <- rbind(c(  0,   0, 0.4, 0.5),
#'               c(  0,   0, 0.3, 0.1),
#'               c(  0,   0,   0,   0),
#'               c(  0,   0,   0,   0))
#'
#' makeLifeTable(matU, startLife = 1, nSteps = 100)
#' makeLifeTable(matU, matF, startLife = 1, nSteps = 100)
#' makeLifeTable(matU, matF, matC, startLife = 1, nSteps = 100)
#' 
#' @export makeLifeTable
makeLifeTable <- function(matU, matF = NULL, matC = NULL, startLife = 1,
                          nSteps = 1000) {
  
  # validate arguments
  if (dim(matU)[1] != dim(matU)[2]) {
    stop("Matrix population model is not a square matrix")
  }
  if (any(colSums(matU) > 1)) {
    warning("matU has at least one stage-specific survival value > 1")
  }
  if (!is.null(matF)) {
    if (any(is.na(matF))) {
      matF[is.na(matF)] = 0
      warning("NAs exist in matF. These have been zero-ed")
    }
  }
  if (!is.null(matC)) {
    if (any(is.na(matC))) {
      matC[is.na(matC)] = 0
      warning("NAs exist in matC. These have been zero-ed")
    }
  }
  
  #Age-specific survivorship (lx)
  lx <- ageSpecificSurv(matU, startLife, nSteps-1)
  
  #Proportion of original cohort dying during each age
  dx = c(lx[1:(length(lx) - 1)] - lx[2:length(lx)], NA)
  
  #Force of mortality
  qx = dx / lx
  
  #Average proportion of individuals alive at the middle of a given age
  Lx = (lx[1:(length(lx) - 1)] + lx[2:length(lx)]) / 2
  Lx[nSteps] <- NA
  
  #Total number of individuals alive at a given age and beyond
  Tx <- sapply(seq_along(Lx), function(x) sum(Lx[x:length(Lx)], na.rm = T))
  Tx[nSteps] <- NA
  
  #Mean life expectancy conditional to a given age
  ex = Tx / lx
  ex[is.infinite(ex)] <- NA
  ex[nSteps] <- NA
  
  #Start to assemble output object
  out = data.frame(
    x = 0:(nSteps - 1),
    lx = lx,
    dx = dx,
    qx = qx,
    Lx = Lx,
    Tx = Tx,
    ex = ex
  )
  
  if (!is.null(matF)) {
    if (sum(matF, na.rm = TRUE) == 0) {
      warning("matF contains only zeros")
    }
    
    #Age-specific fertility (mx)
    out$mx <- ageSpecificRepro(matU, matF, startLife, nSteps-1)
    out$lxmx <- out$lx * out$mx
  }
  
  if (!is.null(matC)) {
    if (sum(matC, na.rm = TRUE) == 0) {
      warning("matC contains only zeros")
    }
    
    #Age-specific clonality (cx)
    out$cx <- ageSpecificRepro(matU, matC, startLife, nSteps-1)
    out$lxcx <- out$lx * out$cx
  }
  
  if (!is.null(matF) & !is.null(matC)) {
    out$mxcx <- out$mx + out$cx
    out$lxmxcx <- out$lx * out$mxcx
  }
  
  return(out)
}
