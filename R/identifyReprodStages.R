#' Identify which stages in a matrix population model are reproductive
#' 
#' Takes a fecundity matrix and returns a vector of logical values (TRUE/FALSE)
#' indicating which stages are reproductive (i.e. exhibit any positive level of
#' reproduction). This is a preparatory step to collapsing the matrix model into
#' a standardized set of stage classes.
#'
#' @param matR The reproductive component of a matrix population model (i.e. a
#'   square projection matrix reflecting transitions due to reproduction; either
#'   sexual, clonal, or both)
#' @param na.handling One of \code{"return.na"}, \code{"return.true"}, or
#'   \code{"return.false"}. Determines how values of \code{NA} within
#'   \code{matR} should be handled. See Value for more details.
#' @return A logical vector of length \code{ncol(matR)}, with values of
#'   \code{FALSE} corresponding to non-reproductive stages and values of
#'   \code{TRUE} corresponding to reproductive stages.\cr\cr For a given matrix
#'   stage (i.e. column of \code{matR}), if there are any positive values of
#'   fecundity, the function will return \code{TRUE}. However, for a given
#'   stage, if there are no positive values of fecundity and one or more values
#'   of \code{NA}, the function will return \code{NA} if \code{na.handling ==
#'   "return.na"}, \code{TRUE} if \code{na.handling == "return.true"}, or
#'   \code{FALSE} if \code{na.handling == "return.false"}.
#' @author Rob Salguero-Gómez <rob.salguero@@zoo.ox.ac.uk>
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' @examples
#' matR1 <- rbind(c( 0, 0.2,   0, 0.5),
#'                c( 0, 0.3,   0, 0.6),
#'                c( 0,   0,   0,   0),
#'                c( 0,   0,   0,   0))
#' 
#' matR2 <- rbind(c(NA,  NA,  NA, 1.1),
#'                c( 0,   0, 0.3, 0.7),
#'                c( 0,   0,   0,   0),
#'                c( 0,   0,   0,   0))
#'
#' identifyReproStages(matR1)
#' 
#' # compare different methods for handling NA
#' identifyReproStages(matR2, na.handling = "return.na")
#' identifyReproStages(matR2, na.handling = "return.true")
#' identifyReproStages(matR2, na.handling = "return.false")
#' @export identifyReproStages
identifyReproStages <- function(matR, na.handling = "return.true") {
  
  # validate arguments
  checkValidMat(matR, fail_any_na = FALSE)
  if (!na.handling %in% c("return.na", "return.true", "return.false")) {
    stop(paste("Argument na.handling must be either 'return.na',",
               "'return.true' or 'return.false'"), call. = FALSE)
  }
  
  if (!any(is.na(matR))) {
    reproStages <- apply(matR, 2, function(x) ifelse(any(x > 0), TRUE, FALSE))
  } else if (na.handling == "return.na") {
    # works because of how function `any` handles NA
    # any(c(0, NA, 0) > 0) will return NA
    # any(c(0, NA, 1) > 0) will return TRUE
    reproStages <- apply(matR, 2, function(x) ifelse(any(x > 0), TRUE, FALSE))
  } else if (na.handling == "return.true") {
    matR[which(is.na(matR))] <- Inf
    reproStages <- apply(matR, 2, function(x) ifelse(any(x > 0), TRUE, FALSE))
  } else if (na.handling == "return.false") {
    matR[which(is.na(matR))] <- 0
    reproStages <- apply(matR, 2, function(x) ifelse(any(x > 0), TRUE, FALSE))
  }
  
  return(reproStages)
}
