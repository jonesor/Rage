#' Calculate generation time from a matrix population model
#'
#' Calculate generation time from a matrix population model, where generation
#' time is defined as the time required for a population to increase by a factor
#' of R0 (the net reproductive rate). See section 5.3.5 of Caswell (2001).
#'
#' @param matU The survival component of a matrix population model (i.e. a
#'   square projection matrix reflecting survival-related transitions; e.g.
#'   progression, stasis, and retrogression)
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' 
#' @details
#' There are multiple definitions of generation time. Here we use \code{log(R0)
#' / log(lambda)}, where \code{R0} is the net reproductive rate (the
#' per-generation population growth rate; Caswell 2001, Sec. 5.3.4), and
#' \code{lambda} is the population growth rate per unit time (the dominant
#' eigenvalue of \code{matU + matR}).
#' 
#' @return Returns generation time. If \code{matU} is singular (often indicating
#'   infinite life expectancy), returns \code{NA}.
#' 
#' @author Patrick Barks <patrick.barks@@gmail.copm>
#' 
#' @references Caswell, H. (2001) Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968
#' 
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
#' # calculate generation time
#' gen_time(matU, matF)
#' 
#' @importFrom popbio lambda
#' @export gen_time
gen_time <- function(matU, matR) {
  
  # leave arg validation to net_repro_rate
  R0 <- net_repro_rate(matU, matR)
  lam <- popbio::lambda(matU + matR)
  
  return(log(R0) / log(lam))
}
