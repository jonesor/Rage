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
#' @author Roberto Salguero-Gómez <rob.salguero@zoo.ox.ac.uk> Hal Caswell
#'   <h.caswell@uva.nl> Owen R. Jones <jones@biology.sdu.dk>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
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
#' matU <- matrix(c(0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0.1), nrow = 4, byrow = T)
#' matF <- matrix(c(0, 0, 5, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.), nrow = 4, byrow = T)
#' matC <- matrix(c(0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4, byrow = T)
#'
#' makeLifeTable(matU, matF, matC, startLife = 1, nSteps = 100)
#'
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
    
    matDim = ncol(matU)
    
    #Age-specific survivorship (lx) (See top function on page 120 in Caswell 2001):
    matUtemp = matU
    survivorship = array(NA, dim = c(nSteps, matDim))
    for (o in 1:nSteps) {
      survivorship[o,] = colSums(matUtemp %*% matU)
      matUtemp = matUtemp %*% matU
    }
    
    lx = survivorship[, startLife]
    lx = c(1, lx[1:(length(lx) - 1)])
    
    #Proportion of original cohort dying during each age
    dx = c(lx[1:(length(lx) - 1)] - lx[2:length(lx)], NA)
    
    #Force of mortality
    qx = dx / lx
    
    #Average proportion of individuals alive at the middle of a given age
    Lx = (lx[1:(length(lx) - 1)] - lx[2:length(lx)]) / 2
    Lx[nSteps] <- NA
    
    #Total number of individuals alive at a given age and beyond
    Tx = cumsum(Lx)
    Tx[nSteps] <- NA
    
    #Mean life expectancy conditional to a given age:
    ex = Tx / lx
    
    for (j in 1:(nSteps - 1)) {
      if (is.infinite(ex[j]))
        ex[j] = NA
    }
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
      if (sum(matF, na.rm = T) == 0) {
        warning("matF contains only 0 values")
      }
      #Age-specific fertility (mx, Caswell 2001, p. 120)
      ageFertility = array(0, dim = c(nSteps, matDim))
      fertMatrix = array(0, dim = c(nSteps, matDim))
      matUtemp2 = matU
      e = matrix(rep(1, matDim))
      for (q in 1:nSteps) {
        fertMatrix = matF %*% matUtemp2 * (as.numeric((solve(
          diag(t(e) %*% matUtemp2)
        ))))
        ageFertility[q,] = colSums(fertMatrix)
        matUtemp2 = matUtemp2 %*% matU
      }
      mx = ageFertility[, startLife]
      mx = c(0, mx[1:(length(mx) - 1)])
      mx[is.nan(mx)] = 0
      out$mx = mx
      out$lxmx = lx * mx
      
      #Net reproductive value
      out$R0Fec = sum(out$lxmx)
      
      #Generation time
      out$TcFec = sum(out$x * out$lxmx) / sum(out$lxmx)
    }
    
    if (!missing(matC)) {
      if (sum(matC, na.rm = T) == 0) {
        warning("matC contains only 0 values")
      }
      #Age-specific clonality (cx)
      ageClonality = array(0, dim = c(nSteps, matDim))
      clonMatrix = array(0, dim = c(nSteps, matDim))
      matUtemp2 = matU
      e = matrix(rep(1, matDim))
      for (q in 1:nSteps) {
        clonMatrix = matC %*% matUtemp2 * (as.numeric((ginv(
          diag(t(e) %*% matUtemp2)
        ))))
        ageClonality[q,] = colSums(clonMatrix)
        matUtemp2 = matUtemp2 %*% matU
      }
      cx = ageClonality[, startLife]
      cx = c(0, cx[1:(length(cx) - 1)])
      cx[is.nan(cx)] = 0
      out$cx = cx
      out$lxcx = lx * cx
      
      #Net reproductive value
      out$R0Clo = sum(out$lxcx)
      
      #Generation time
      out$TcClo = sum(out$x * out$lxcx) / sum(out$lxcx)
    }
    
    if (!missing(matF) & !missing(matC)) {
      mxcx = mx + cx
      out$mxcx = mxcx
      out$lxmxcx = lx * mxcx
      
      #Net reproductive value
      out$R0FecClo = sum(out$lxmxcx)
      
      #Generation time
      out$TcFecClo = sum(out$x * out$lxmxcx) / sum(out$lxmxcx)
    }
    
    return(as.life.table(out))
  }
