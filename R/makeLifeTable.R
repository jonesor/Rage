#' Produces a life table from a matrix population model
#'
#' This function uses age-from-stage decompositions to calculate a life table
#' from a matrix population model
#'
#' A detailed description of these methods can be found in section 5.3 of
#' Caswell (2001) and the supplementary information of Jones et al. (2014).
#'
#' @param matU A matrix containing only survival-dependent processes (e.g.
#'   progression, stasis, retrogression).
#' @param matF A matrix containing only sexual reproduction, with zeros
#'   elsewhere.
#' @param matC A matrix containing only clonal reproduction, with zeros
#'   elsewhere. If not provided, it defaults to a matrix with all zeros.
#' @param startLife The first stage at which the author considers the beginning
#'   of life in the life cycle of the species. It defaults to the first stage.
#' @param nSteps Number of time steps for which the life table will be
#'   calculated. This is on the same units as the matrix population model - see
#'   MatrixPeriodicity in metadata of COMPADRE/COMADRE.
#' @return A `data.frame` with between 5 and 8 columns. If just `matU` is
#'   provided, the `data.frame` has columns `x` (age), `lx` (survivorship), 'dx'
#'   (proportion of the original cohort dying during each stage), 'qx' (force of
#'   mortality), and 'ex' (age-specific life expectancy). If `matF` is provided
#'   (in addition to `matU`), `mx` (age-specific sexual reproduction) and lxmx
#'   (# sexually-produced recruits per capita in each age) are included in the
#'   output. If `matC` is provided (in addition to `matU`), `cx` (age-specific
#'   clonal reproduction) and lxcx (# clonally-produced recruits per capita in
#'   each age) are included in the output.
#' @note %% ~~further notes~~
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk> 
#' @author Hal Caswell <h.caswell@@uva.nl> 
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @references
#'
#' Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and
#' Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968
#'
#' Caswell, H. (2006) Applications of Markov chains in demography. pp. 319-334
#' in A.N. Langville and W.J. Stewart (editors) MAM2006: Markov Anniversary
#' Meeting. Boson Books, Raleigh, North Caroline, USA
#'
#' Jones, O.R. et al. (2014) Diversity of ageing across the tree of life.
#' Nature, 505(7482), 169–173
#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#'
#' matU <- matrix(c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = TRUE)
#' matF <- matrix(c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.), nrow = 4, byrow = TRUE)
#' matC <- matrix(c(0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = TRUE)
#'
#' makeLifeTable(matU, matF, matC, startLife = 1, nSteps = 100)
#' @export makeLifeTable
makeLifeTable <-
  function(matU,
           matF = NULL,
           matC = NULL,
           startLife = 1,
           nSteps = 1000) {
    #Error checks
    if (dim(matU)[1] != dim(matU)[2])
      stop("Your matrix population model is not a square matrix")
    if (length(which(colSums(matU) > 1)) > 0)
      warning("matU has at least one stage-specific survival value > 1")
    
    if (!missing(matF)) {
      if (any(is.na(matF))) {
        matF[is.na(matF)] = 0
        warning("NAs exist in matF. These have been zero-ed")
      }
    }
    if (!missing(matC)) {
      if (any(is.na(matC))) {
        matC[is.na(matC)] = 0
        warning("NAs exist in matC. These have been zero-ed")
      }
    }
    
    #Age-specific survivorship (lx):
    lx <- ageSpecificSurv(matU, startLife, nSteps-1)
    # lx <- head(lx, -1) # remove lx[nSteps+1]
    
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
    
    #Mean life expectancy conditional to a given age:
    ex = Tx / lx
    ex[is.infinite(ex)] <- NA
    ex[nSteps] <- NA
    
    #Start to assemble output object
    out = list(
      x = 0:(nSteps - 1),
      lx = lx,
      dx = dx,
      qx = qx,
      Lx = Lx,
      Tx = Tx,
      ex = ex
    )
    
    if (!missing(matF)) {
      if (sum(matF, na.rm = TRUE) == 0) {
        warning("matF contains only 0 values")
      }
      
      #Age-specific fertility (mx)
      out$mx <- ageSpecificRepro(matU, matF, startLife, nSteps-1)
      # out$mx <- head(out$mx, -1) # remove mx[nSteps+1]
      out$lxmx <- out$lx * out$mx
      
      #Net reproductive value
      out$R0Fec = sum(out$lxmx)
      
      #Generation time
      out$TcFec = sum(out$x * out$lxmx) / sum(out$lxmx)
    }
    
    if (!missing(matC)) {
      if (sum(matC, na.rm = TRUE) == 0) {
        warning("matC contains only 0 values")
      }
      
      #Age-specific clonality (cx)
      out$cx <- ageSpecificRepro(matU, matC, startLife, nSteps-1)
      # out$cx <- head(out$cx, -1) # remove cx[nSteps+1]
      out$lxcx <- out$lx * out$cx
      
      #Net reproductive value
      out$R0Clo = sum(out$lxcx)
      
      #Generation time
      out$TcClo = sum(out$x * out$lxcx) / sum(out$lxcx)
    }
    
    if (!missing(matF) & !missing(matC)) {
      out$mxcx = out$mx + out$cx
      out$lxmxcx = out$lx * out$mxcx
      
      #Net reproductive value
      out$R0FecClo = sum(out$lxmxcx)
      
      #Generation time
      out$TcFecClo = sum(out$x * out$lxmxcx) / sum(out$lxmxcx)
    }
    
    return(out)
  }
