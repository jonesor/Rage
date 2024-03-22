#' Calculate Keyfitz entropy for a matrix population model
#'
#' Computes Keyfitz entropy from the U submatrix of a matrix population model.
#'
#' @param m A square numeric matrix representing the U submatrix of a matrix
#'   population model. For age-based matrices, the dimension corresponds to the
#'   number of age classes
#' @param type An argument specifying the type of matrix model used in the
#'   calculations. This is necessary because the calculations for age vs.
#'   stage-based matrices are different. Possible options are `age` and `stage`.
#'   The latter is not yet implemented. Defaults to `age`.
#' 
#' @return Returns a single numeric value representing the Keyfitz entropy
#' for the given matrix. This value quantifies the dispersion of age at death.
#'
#' @author Lotte de Vries <c.devries@@uva.nl>
#'
#' @references Keyfitz, N. 1977. Applied Mathematical Demography. New York:
#'   Wiley.
#'   
#'   de Vries, C., Bernard, C., & Salguero-Gómez, R. 2023. Discretising
#'   Keyfitz' entropy for studies of actuarial senescence and comparative
#'   demography. Methods in Ecology and Evolution, 14, 1312–1319.
#'   https://doi.org/10.1111/2041-210X.14083
#'   
#' @family life history traits
#' @examples
#' data(leslie_mpm1)
#' 
#' entropy_k2(leslie_mpm1$matU)
#'
#' @export

entropy_k2 <- function(m, type = "age") {
  if (!inherits(m, "matrix")) {
    stop("m must be a matrix")
  }
  if (!is.numeric(m)) {
    stop("m must be a numeric matrix")
  }
  if (nrow(m) != ncol(m)) {
    stop("m must be a square matrix")
  }
if(type == "age"){
  maxage <- dim(m)[1]
  e0 <- rep(0, maxage)
  e0[1] <- 1
  l <- rep(0, 2 * maxage)
  temp <- rep(1, maxage) - colSums(m)
  M <- diag(temp)
  N <- solve(diag(maxage) - m)
  eta1 <- rep(1, maxage) %*% N
  B <- M %*% N
  etadagger <- eta1 %*% B
  eta0 <- eta1 %*% e0
  H_new <- etadagger[1] / eta0
  
  return(as.vector(H_new))
}
  if(type == "stage"){
    stop("This has not yet been implemented.")
  }
}
