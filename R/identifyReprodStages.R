#' Identify which stages in a matrix model are reproductive
#' 
#' Takes a fecundity matrix and returns a vector of logical values (TRUE/FALSE)
#' indicating which stages are reproductive (i.e. exhibit any positive level of
#' reproduction). This is a preparatory step to collapsing the matrix model into
#' combined propagule, pre-reproductive, reproductive and post-reproductive
#' stages.
#'
#' @export
#' @param matF A fecundity matrix
#' @param na.handling One of \code{"return.na"}, \code{"return.true"}, or
#'   \code{"return.false"}. Determines how values of \code{NA} within
#'   \code{matF} should be handled. See Value for more details.
#' @return A logical vector of length \code{ncol(matF)}, with values of
#'   \code{FALSE} corresponding to non-reproductive stages and values of
#'   \code{TRUE} corresponding to reproductive stages.\cr\cr For a given matrix
#'   stage (i.e. column of \code{matF}), if there are any positive values of
#'   fecundity, the function will return \code{TRUE}. However, for a given
#'   stage, if there are no positive values of fecundity and one or more values
#'   of \code{NA}, the function will return \code{NA} if \code{na.handling ==
#'   "return.na"}, \code{TRUE} if \code{na.handling == "return.true"}, or
#'   \code{FALSE} if \code{na.handling == "return.false"}.
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' @examples
#' matF1 <- rbind(c(0, 0.2, 0, 0.5, 0), c(0, 0.3, 0.2, 0.6, 0), c(0, 0, 0, 0,
#' 0), c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
#'
#' matF2 <- rbind(c(NA, NA, NA, 1.11), c(0, 0, 0.31, 0.78), c(0, 0, 0, 0), c(0,
#' 0, 0, 0), c(0, 0, 0, 0))
#'
#' identifyReproStages(matF1)
#' 
#' # compare different methods for handling NA
#' identifyReproStages(matF2, na.handling = "return.na")
#' identifyReproStages(matF2, na.handling = "return.true")
#' identifyReproStages(matF2, na.handling = "return.false")
#' 
#' \dontrun{
#' # invalid setting for argument na.handling
#' identifyReproStages(matF2, na.handling = "NA")
#' }
identifyReproStages <- function(matF, na.handling = "return.true") {
  if (!na.handling %in% c("return.na", "return.true", "return.false")) {
    stop("Argument na.handling must be either 'return.na', 'return.true', or 'return.false'")
  } else if (all(is.na(matF))) {
    stop("All elements of matF are NA")
  } else if (!any(is.na(matF))) {
    reproStages <- apply(matF, 2, function(x) ifelse(any(x > 0), TRUE, FALSE))
  } else if (na.handling == "return.na") {
    # works because of how function `any` handles NA
    # any(c(0, NA, 0) > 0) will return NA
    # any(c(0, NA, 1) > 0) will return TRUE
    reproStages <- apply(matF, 2, function(x) ifelse(any(x > 0), TRUE, FALSE))
  } else if (na.handling == "return.true") {
    matF[which(is.na(matF))] <- Inf
    reproStages <- apply(matF, 2, function(x) ifelse(any(x > 0), TRUE, FALSE))
  } else if (na.handling == "return.false") {
    matF[which(is.na(matF))] <- 0
    reproStages <- apply(matF, 2, function(x) ifelse(any(x > 0), TRUE, FALSE))
  }
  
  return(reproStages)
}
