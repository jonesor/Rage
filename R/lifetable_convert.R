#' Convert between age-specific survivorship, survival, or mortality hazard
#'
#' Convert between vectors of age-specific survivorship (lx), survival
#' probability (px), or mortality hazard (hx). Input vectors must be arranged in
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
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @family {life tables}
#' 
#' @references Ergon, T., Borgan, Ø., Nater, C. R., & Vindenes, Y. 2018. The
#'   utility of mortality hazard rates in population analyses. Methods in
#'   Ecology and Evolution, 9, 2046-2056. <doi:10.1111/2041-210X.13059>
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
