#' Convert between age-specific survivorship, survival, or mortality hazard
#'
#' Convert between vectors of age-specific survivorship (\code{lx}), survival
#' probability (\code{px}), or mortality hazard (\code{hx}). Input vectors must be arranged in
#' order of increasing age, starting with age 0.
#'
#' @param lx Vector of age-specific survivorship.
#' @param px Vector of age-specific survival probabilities.
#' @param hx Vector of age-specific mortality hazards.
#' 
#' @return A vector.
#' 
#' @details \code{lx} gives the proportional survivorship to the start of age
#'   class \code{x} (where survivorship at first age class is defined as 1),
#'   \code{px} gives the probability of survival between age \code{x} and
#'   \code{x+1}, and \code{hx} gives the time-averaged mortality hazard (also
#'   called force of mortality) between age \code{x} and \code{x+1}.
#'   
#' @note Note that the units of time for the returned vectors (i.e., `x`) are
#'   the same as the projection interval (`ProjectionInterval`) of the MPM.
#'   
#' @author Patrick Barks <patrick.barks@@gmail.com>
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
#' Ergon, T., Borgan, Ø., Nater, C. R., & Vindenes, Y. 2018. The utility of
#' mortality hazard rates in population analyses. Methods in Ecology and
#' Evolution, 9, 2046-2056. <doi:10.1111/2041-210X.13059>
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
#' lx <- c(1, 0.8, 0.7, 0.5, 0.3, 0.1)
#' 
#' # convert from lx
#' px <- lx_to_px(lx)
#' hx <- lx_to_hx(lx)
#' 
#' # convert from px
#' lx <- px_to_lx(px)
#' hx <- px_to_hx(px)
#' 
#' # convert from hx
#' lx <- hx_to_lx(hx)
#' px <- hx_to_px(hx)
#' 
#' @name lifetable_convert
NULL



#' @rdname lifetable_convert
#' @export lx_to_px
lx_to_px <- function(lx) {
  if (length(lx) == 1) {
    NA_real_
  } else {
    c(lx[-1] / lx[-length(lx)], NA)
  }
}


#' @rdname lifetable_convert
#' @export lx_to_hx
lx_to_hx <- function(lx) {
  if (length(lx) == 1) {
    NA_real_
  } else {
    px <- lx_to_px(lx)
    px_to_hx(px)
  }
}


#' @rdname lifetable_convert
#' @export px_to_lx
px_to_lx <- function(px) {
  c(1, cumprod(px[-length(px)]))
}


#' @rdname lifetable_convert
#' @export px_to_hx
px_to_hx <- function(px) {
  -log(px)
}


#' @rdname lifetable_convert
#' @export hx_to_lx
hx_to_lx <- function(hx) {
  px <- hx_to_px(hx)
  px_to_lx(px)
}


#' @rdname lifetable_convert
#' @export hx_to_px
hx_to_px <- function(hx) {
  exp(-hx)
}
