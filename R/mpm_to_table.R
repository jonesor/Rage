#' Generate a life table from a matrix population model
#'
#' This function uses age-from-stage decomposition methods to generate a life
#' table from a matrix population model. A detailed description of these methods
#' can be found in section 5.3 of Caswell (2001).
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matF (Optional) The sexual component of a matrix population model
#'   (i.e. a square projection matrix reflecting transitions due to sexual
#'   reproduction)
#' @param matC (Optional) The clonal component of a matrix population model
#'   (i.e. a square projection matrix reflecting transitions due to clonal
#'   reproduction).
#' @param start The index of the first stage at which the author considers the
#'   beginning of life. Defaults to 1. Alternately, a numeric vector giving the
#'   starting population vector (in which case \code{length(start)} must match
#'   \code{ncol(matU))}. See section \emph{Starting from multiple stages}.
#' @param xmax Maximum age to which the life table will be calculated (defaults
#'   to \code{1000}). Time steps are in the same units as the matrix population
#'   model (see MatrixPeriodicity column in metadata of COM(P)ADRE).
#' @param lx_crit Minimum value of lx to which age-specific traits will be
#'   calculated (defaults to \code{1e-4}).
#' 
#' @return A \code{data.frame} containing 7-13 columns.
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
#' 
#' @section Starting from multiple stages:
#' Rather than specifying argument \code{start} as a single stage class from
#' which all individuals start life, it may sometimes be desirable to allow for
#' multiple starting stage classes. For example, if we want to start our
#' calculation of age-specific traits from reproductive maturity (i.e. first
#' reproduction), we should account for the possibility that there may be
#' multiple stage classes in which an individual could first reproduce.
#' 
#' To specify multiple starting stage classes, specify argument \code{start} as
#' the desired starting population vector (\strong{n1}), giving the proportion
#' of individuals starting in each stage class (the length of \code{start}
#' should match the number of columns in the relevant MPM).
#' 
#' See function \code{\link{mature_distrib}} for calculating the proportion of
#' individuals achieving reproductive maturity in each stage class.
#' 
#' @note The life table is calculated recursively until the age class (x)
#'   reaches \code{xmax} or survivorship (lx) falls below \code{lx_crit} —
#'   whichever comes first. To force calculation to \code{xmax}, set
#'   \code{lx_crit = 0}. Conversely, to force calculation to \code{lx_crit}, set
#'   \code{xmax = Inf}.
#' 
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk> 
#' @author Hal Caswell <h.caswell@@uva.nl> 
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' 
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
#' 
#' @examples
#' data(mpm1)
#' 
#' mpm_to_table(matU = mpm1$matU, start = 2, xmax = 15)
#' mpm_to_table(matU = mpm1$matU, matF = mpm1$matF, start = 2, xmax = 15)
#' 
#' ### starting from first reproduction
#' repstages <- repro_stages(mpm1$matF)
#' n1 <- mature_distrib(matU = mpm1$matU, start = 2, repro_stages = repstages)
#' mpm_to_table(matU = mpm1$matU, start = n1)
#' 
#' @export mpm_to_table
mpm_to_table <- function(matU, matF = NULL, matC = NULL, start = 1L,
                         xmax = 1000, lx_crit = 1e-4) {
  
  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  if (!is.null(matF)) checkValidMat(matF)
  if (!is.null(matC)) checkValidMat(matC)
  checkValidStartLife(start, matU, start_vec = TRUE)
  
  #Age-specific survivorship (lx)
  lx <- mpm_to_lx(matU, start, xmax, lx_crit)
  
  N <- length(lx)
  
  #Proportion of original cohort dying during each age
  dx <- c(lx[1:(N - 1)] - lx[2:N], NA)
  
  #Force of mortality
  qx <- dx / lx
  
  #Average proportion of individuals alive at the middle of a given age
  Lx <- (lx[1:(N - 1)] + lx[2:N]) / 2
  Lx[N] <- NA
  
  #Total number of individuals alive at a given age and beyond
  TxFn <- function(x) sum(Lx[x:length(Lx)], na.rm = TRUE)
  Tx <- vapply(seq_along(Lx), TxFn, numeric(1))
  Tx[N] <- NA
  
  #Mean life expectancy conditional to a given age
  ex <- Tx / lx
  ex[is.infinite(ex)] <- NA
  ex[N] <- NA
  
  #Start to assemble output object
  out <- data.frame(
    x = 0:(N - 1),
    lx = lx,
    dx = dx,
    qx = qx,
    Lx = Lx,
    Tx = Tx,
    ex = ex
  )
  
  if (!is.null(matF)) {
    out$mx <- mpm_to_mx(matU, matF, start, N-1)
    out$lxmx <- out$lx * out$mx
  }
  
  if (!is.null(matC)) {
    out$cx <- mpm_to_mx(matU, matC, start, N-1)
    out$lxcx <- out$lx * out$cx
  }
  
  if (!is.null(matF) & !is.null(matC)) {
    out$mxcx <- out$mx + out$cx
    out$lxmxcx <- out$lx * out$mxcx
  }
  
  return(out)
}
