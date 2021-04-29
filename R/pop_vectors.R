#' Derive a hypothetical set of population vectors corresponding to a
#' time-series of matrix population models
#'
#' Derive a hypothetical set of population vectors (i.e. population size
#' distributions across stages) given a time-series of matrix population models
#' (MPMs), by taking the stable stage distribution of the mean matrix as the starting
#' vector (or optionally, a uniform or random starting vector), and deriving
#' subsequent vectors through recursive population projection.
#' 
#' This function is useful for providing population vectors as input to the
#' \code{\link{perturb_stochastic}} function which calculates stochastic
#' elasticities given a time-series of matrix population models and
#' corresponding population vectors, using the method described in Haridas et
#' al. (2009).
#' 
#' @param A A list of MPMs (i.e. square population projection matrices).
#' @param start Method to derive the first population vector in the series.
#'   Either `stable.stage` to use the stable stage distribution of the mean
#'   matrix as the starting vector, `uniform` to use a uniform starting vector
#'   (all elements equal), or `random` to use a randomly-generated starting
#'   vector. Defaults to `stable.stage`.
#' 
#' @return A list of population vectors
#' 
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @family perturbation analysis
#' 
#' @references Haridas, C. V., Tuljapurkar, S., & Coulson, T. 2009. Estimating
#'   stochastic elasticities directly from longitudinal data. Ecology Letters,
#'   12, 806-812. <doi:10.1111/j.1461-0248.2009.01330.x>
#' 
#' @examples 
#' # generate list of matrices
#' matA_l <- replicate(5, matrix(runif(9), 3, 3), simplify = FALSE)
#' 
#' # calculate corresponding population vectors
#' pop_vectors(matA_l)
#' pop_vectors(matA_l, start = "uniform")
#' pop_vectors(matA_l, start = "random")
#' 
#' @importFrom popbio stable.stage
#' @importFrom stats runif
#' @export pop_vectors
pop_vectors <- function(A, start = "stable.stage") {
  
  A_mean <- meanMat(A)
  s <- nrow(A_mean)
  w <- list()
  
  if (start == "stable.stage") {
    w[[1]] <- stable.stage(A_mean)
  } else if (start == "uniform") {
    w[[1]] <- rep(1/s, s)
  } else if (start == "random") {
    w[[1]] <- runif(s)
    w[[1]] <- w[[1]] / sum(w[[1]])
  } else {
    stop("Argument start must be one of 'stable.stage', 'uniform', or 'random'",
         .call = FALSE)
  }
  
  for (i in 2:length(A)) {
    w[[i]] <- A[[i-1]] %*% w[[i-1]]
    w[[i]] <- as.numeric(w[[i]] / sum(w[[i]]))
  }
  
  return(w)
}
