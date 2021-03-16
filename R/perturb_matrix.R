#' Perturbation analysis of a matrix population model
#'
#' Perturbs elements within a matrix population model and measures the response 
#' (sensitivity or elasticity) of the per-capita population growth rate at 
#' equilibrium (\eqn{\lambda}), or, with a user-supplied function, any other 
#' demographic statistic.
#'
#' @param matA A matrix population model (i.e. a square projection matrix).
#' @param pert Magnitude of the perturbation. Defaults to \code{1e-6}.
#' @param type Whether to return "sensitivity" or "elasticity" values.
#' @param demog_stat The demographic statistic to be used, as in "the
#'   sensitivity/elasticity of \code{demog_stat} to matrix element 
#'   perturbations." Defaults to the per-capita population growth rate at 
#'   equilibrium (\eqn{\lambda}). Also accepts a user-supplied function that 
#'   performs a calculation on a projection matrix and returns a single numeric 
#'   value.
#' @param ... Additional arguments passed to the function \code{demog_stat}
#' 
#' @return A sensitivity or elasticity matrix.
#' 
#' @author Rob Salguero-Gomez <rob.salguero@@zoo.ox.ac.uk>
#' 
#' @family {perturbation analysis}
#' 
#' @references Caswell, H. 2001. Matrix Population Models: Construction,
#'   Analysis, and Interpretation. Sinauer Associates; 2nd edition. ISBN:
#'   978-0878930968 
#'   
#' @examples
#' matA <- rbind(c(0.1,   0, 1.5, 4.6),
#'               c(0.5, 0.2, 0.1,   0),
#'               c(  0, 0.3, 0.3, 0.1),
#'               c(  0,   0, 0.5, 0.6))
#'
#' perturb_matrix(matA)
#' 
#' # use a larger perturbation than the default
#' perturb_matrix(matA, pert = 0.01)
#' 
#' # calculate the sensitivity/elasticity of the damping ratio to perturbations
#' damping <- function(matA) {  # define function for damping ratio
#'   eig <- eigen(matA)$values
#'   dm <- rle(Mod(eig))$values
#'   return(dm[1] / dm[2])
#' }
#' 
#' perturb_matrix(matA, demog_stat = "damping")
#' 
#' @importFrom popbio lambda
#' @export perturb_matrix
perturb_matrix <- function(matA, pert = 1e-6, type = "sensitivity",
                           demog_stat = "lambda", ...) {
  
  # validate arguments
  checkValidMat(matA)
  
  # get statfun
  if (is.character(demog_stat) && demog_stat == "lambda") {
    statfun <- popbio::lambda
  } else {
    statfun <- try(match.fun(demog_stat), silent = TRUE)
    if (class(statfun) == "try-error") {
      stop("demog_stat must be 'lambda' or the name of a function that ",
           "returns a single numeric value", call. = FALSE)
    }
  }

  # matrix dimension
  m <- nrow(matA)
  
  # statfun
  stat <- statfun(matA, ...)
  
  # initialize sensitivity matrix
  sensA <- matrix(NA, ncol = m, nrow = m)
  
  # matrix perturbation
  for (i in 1:m) {
    for (j in 1:m) {
      pertA <- matA
      pertA[i, j] <- pertA[i, j] + pert
      statPert <- statfun(pertA, ...)
      sensA[i, j] <- (stat - statPert) / (matA[i, j] - pertA[i, j])
    }
  }
  
  # keep only real values
  sensA <- Re(sensA)
  
  # copy attributes from matA
  attributes(sensA) <- attributes(matA)
  
  # return
  if (type == "sensitivity") {
    return(sensA)
  } else {
    return(sensA * matA / stat)
  }
}

