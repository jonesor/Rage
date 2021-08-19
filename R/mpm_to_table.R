#' Generate a life table from a matrix population model
#'
#' This function uses age-from-stage decomposition methods to generate a life
#' table from a matrix population model. A detailed description of these methods
#' can be found in section 5.3 "Age-specific traits from stage-specific models"
#' of Caswell (2001).
#' 
#' @param matU The survival component of a matrix population model (i.e., a
#'   square projection matrix reflecting survival-related transitions; e.g.,
#'   progression, stasis, and/or retrogression). Optionally with named rows and
#'   columns indicating the corresponding life stage names.
#' @param matF (Optional) The sexual component of a matrix population model
#'   (i.e., a square projection matrix reflecting transitions due to sexual
#'   reproduction). Optionally with named rows and columns indicating the
#'   corresponding life stage names.
#' @param matC (Optional) The clonal component of a matrix population model
#'   (i.e., a square projection matrix reflecting transitions due to clonal
#'   reproduction). Optionally with named rows and columns indicating the
#'   corresponding life stage names.
#' @param start The index (or stage name) of the first stage at which the author
#'   considers the beginning of life. Defaults to \code{1}. Alternately, a numeric
#'   vector giving the starting population vector (in which case
#'   \code{length(start)} must match \code{ncol(matU))}. See section
#'   \emph{Starting from multiple stages}.
#' @param xmax Maximum age to which the life table will be calculated (defaults
#'   to \code{1000}). Time steps are in the same units as the matrix population
#'   model (see MatrixPeriodicity metadata variable COM(P)ADRE).
#' @param lx_crit Minimum value of lx to which age-specific traits will be
#'   calculated (defaults to \code{1e-4}).
#' @param radix The starting number of individuals in the synthetic life table
#'   (defaults to \code{1}). If \code{radix} is set to 1, a simplified life
#'   table is produced.
#'
#' @return A \code{data.frame} containing a variable number columns, depending
#'   on input variables. Columns include:
#'   
#'   \item{x}{age at the start of the age interval \code{[x, x+1)}}
#'   \item{Nx}{The number of individuals alive at age x. The initial number is
#'   set with \code{radix}}
#'   \item{Dx}{proportion of original cohort dying during the age interval \code{[x, x+1)}}
#'   \item{lx}{survivorship, defined as the proportion of initial cohort surviving to the start of
#'   age interval \code{[x, x+1)}}
#'   \item{dx}{proportion of original cohort dying in the age interval \code{[x, x+1)}}
#'   \item{ax}{The average time survived within the interval by those that die
#'   during the age interval \code{[x, x+1)}. Assumed to be 0.5}
#'   \item{hx}{force of mortality (hazard) during the age interval \code{[x, x+1)}}
#'   \item{qx}{probability of death during the interval \code{[x, x+1)} for
#'   those entering the interval}
#'   \item{px}{probability of survival for the interval \code{[x, x+1)} for
#'   those entering the interval}
#'   \item{Lx}{total person-years lived during the interval \code{[x, x+1)}} 
#'   \item{Tx}{total person years lived beyond age x}
#'   \item{ex}{remaining life expectancy at age x}
#'
#' If \code{matF} is provided, also includes:
#'   \item{mx}{per-capita rate of sexual reproduction during the interval
#'    \code{[x, x+1)} }
#'   \item{lxmx}{expected number of sexual offspring per original
#'   cohort member produced during the interval \code{[x, x+1)}}
#'
#' If \code{matC} is provided, also includes:
#'   \item{cx}{per-capita rate of clonal reproduction  during the interval
#'   \code{[x, x+1)}}
#'   \item{lxcx}{expected number of clonal offspring per original
#'   cohort member produced during the interval \code{[x, x+1)}}
#'
#' If both \code{matF} and \code{matC} are provided, also includes:
#'   \item{mxcx}{per-capita rate of total reproduction (sexual + clonal) during
#'   the interval \code{[x, x+1)}}
#'   \item{lxmxcx}{expected number of total offspring (sexual + clonal) per
#'   original cohort member produced during the interval \code{[x, x+1)}}
#'
#' @section Starting from multiple stages: Rather than specifying argument
#'   \code{start} as a single stage class from which all individuals start life,
#'   it may sometimes be desirable to allow for multiple starting stage classes.
#'   For example, if the user wants to start the calculation of age-specific
#'   traits from reproductive maturity (i.e., first reproduction), the user
#'   should account for the possibility that there may be multiple stage classes
#'   in which an individual could first reproduce.
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
#' @note The life table calculations assume that the final age interval is
#'   closed and that all remaining individuals die in this interval. Therefore,
#'   for this interval, the probability of death \code{qx} is 1, the probability
#'   of survival \code{px} is 0 and, because we assume that deaths are evenly
#'   distributed during the interval, the remaining life expectancy for
#'   individuals at the start of the interval is 0.5.
#'   
#'   If \code{lx_crit} is sufficiently small that only a very small proportion
#'   of the cohort reach this age (i.e., < 0.05), this should have minimal impact
#'   on results. Nevertheless, for many analyses, the final row of the life
#'   table should be treated with caution and perhaps removed from subsequent
#'   analyses.
#'   
#' @note Note that the units of time (e.g.. `x` and `ex`) in the returned life
#'   table are the same as the projection interval (`ProjectionInterval`) of the
#'   MPM.
#'   
#' @author Owen R. Jones <jones@@biology.sdu.dk>
#' @author Roberto Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Hal Caswell <h.caswell@@uva.nl>
#'
#' @family life tables
#'
#' @references
#' Caswell, H. 2001. Matrix Population Models: Construction, Analysis, and
#' Interpretation. Sinauer Associates; 2nd edition. ISBN: 978-0878930968
#'
#' Caswell, H. 2006. Applications of Markov chains in demography. pp. 319-334 in
#' A.N. Langville and W.J. Stewart (editors) MAM2006: Markov Anniversary
#' Meeting. Boson Books, Raleigh, North Caroline, USA
#'
#' Horvitz, C. & Tuljapurkar, S. 2008. Stage dynamics, period survival, and
#' mortality plateaus. The American Naturalist 172: 203-2015.
#' <doi:10.1086/589453>
#'
#' Jones, O. R., Scheuerlein, A., Salguero-Gomez, R., Camarda, C. G., Schaible,
#' R., Casper, B. B., Dahlgren, J. P., Ehrlén, J., García, M. B., Menges, E.,
#' Quintana-Ascencio, P. F., Caswell, H., Baudisch, A. & Vaupel, J. 2014.
#' Diversity of ageing across the tree of life. Nature 505, 169-173.
#' <doi:10.1038/nature12789>
#' 
#' Jones O. R. 2021. Life tables: Construction and interpretation In:
#' Demographic Methods Across the Tree of Life. Edited by Salguero-Gomez R &
#' Gamelon M. Oxford University Press. Oxford, UK. ISBN: 9780198838609
#' 
#' Preston, S., Heuveline, P., & Guillot, M. 2000. Demography: Measuring and
#' Modeling Population Processes. Wiley. ISBN: 9781557864512
#'
#' @examples
#' data(mpm1)
#'
#' mpm_to_table(matU = mpm1$matU, start = 2, xmax = 15)
#' mpm_to_table(matU = mpm1$matU, start = "small", xmax = 15) # equivalent using named life stages
#' mpm_to_table(matU = mpm1$matU, matF = mpm1$matF, start = 2, xmax = 15)
#'
#' ### starting from first reproduction
#' repStages <- repro_stages(mpm1$matF)
#' n1 <- mature_distrib(matU = mpm1$matU, start = 2, repro_stages = repStages)
#' mpm_to_table(matU = mpm1$matU, start = n1)
#' @export mpm_to_table
mpm_to_table <- function(matU, matF = NULL, matC = NULL, start = 1L,
                         xmax = 1000, lx_crit = 1e-4, radix = 1) {

  # validate arguments
  checkValidMat(matU, warn_surv_issue = TRUE)
  if (!is.null(matF)) checkValidMat(matF)
  if (!is.null(matC)) checkValidMat(matC)
  checkValidStartLife(start, matU, start_vec = TRUE)

  # Age-specific survivorship (lx)
  lx <- mpm_to_lx(matU, start, xmax, lx_crit)
  if (lx[length(lx)] > 0.05) {
    warning(strwrap(prefix = " ", initial = "", "There are still a large proportion of the synthetic
  cohort remaining alive in the final row of the life table. Consider changing values for `lx_crit` or `xmax`"))
  }
  # Number left alive at age x
  Nx <- lx * radix

  # Number of age groups
  N <- length(lx)

  # Proportion of ORIGINAL cohort dying during each age interval
  dx <- c(lx[1:(N - 1)] - lx[2:N], lx[length(lx)])

  # Number of deaths occuring during the interval (assuming an initial
  # population size defined by the radix)
  # If radix = 1 (the default) then dx = Dx
  Dx <- dx * radix

  # Probability of dying between x and x+1 (qx)
  qx <- dx / lx

  # Probability of surviving between x and x+1 (px)
  px <- 1 - qx


  # Person-years lived between x and x+1
  # We use an ax value of 0.5, which means that we implicitly assume that deaths
  # are distributed uniformly across the interval. ax is defined as the average
  # number of years lived in the interval by those dying in the interval
  
  ax <- 0.5
  Lx <- (px * lx) + (ax * dx)

  # Person=years lived above age x
  Tx <- rev(cumsum(rev(Lx)))

  # Death rate (hazard) in the cohort between ages x and x + 1
  hx <- dx / Lx

  # Life expectancy from age x
  ex <- Tx / lx

  # Force of mortality (hazard)
  qx <- dx / lx

  # Start to assemble output object
  # If radix is set to be a value other than 1 a full life table is produced,
  # otherwise a streamlined version is produced. In most use cases the
  # streamlined version is sufficient.
  if (radix != 1) {
    out <- data.frame(
      x = 0:(N - 1),
      Nx = Nx,
      Dx = Dx,
      lx = lx,
      dx = dx,
      ax = ax,
      hx = hx,
      qx = qx,
      px = px,
      Lx = Lx,
      Tx = Tx,
      ex = ex
    )
  } else {
    out <- data.frame(
      x = 0:(N - 1),
      lx = lx,
      dx = dx,
      hx = hx,
      qx = qx,
      px = px,
      ex = ex
    )
  }

  if (!is.null(matF)) {
    out$mx <- mpm_to_mx(matU, matF, start, N - 1)
    out$lxmx <- out$lx * out$mx
  }

  if (!is.null(matC)) {
    out$cx <- mpm_to_mx(matU, matC, start, N - 1)
    out$lxcx <- out$lx * out$cx
  }

  if (!is.null(matF) & !is.null(matC)) {
    out$mxcx <- out$mx + out$cx
    out$lxmxcx <- out$lx * out$mxcx
  }

  return(out)
}
