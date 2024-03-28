#' Calculate Keyfitz entropy for an age-based matrix population model
#'
#' Computes Keyfitz entropy from the U submatrix of an age-based matrix
#' population model.
#'
#' @param Umat A square numeric matrix representing the U submatrix of a matrix
#'   population model. The dimension corresponds to the number of age classes
#'
#' @return Returns a single numeric value representing the Keyfitz entropy
#' for the given matrix. This value quantifies the dispersion of age at death.
#'
#' @author Lotte de Vries <c.devries@@uva.nl>
#' @author Owen Jones <jones@@biology.sdu.dk>
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
#' entropy_k_age(leslie_mpm1$matU)
#' @export

entropy_k_age <- function(Umat) {
  if (!inherits(Umat, "matrix")) {
    stop("m must be a matrix")
  }

  if (!is.numeric(Umat)) {
    stop("m must be a numeric matrix")
  }
  if (nrow(Umat) != ncol(Umat)) {
    stop("m must be a square matrix")
  }

  # Function from Lotte de Vries
  maxage <- dim(Umat)[1]
  e0 <- rep(0, maxage)
  e0[1] <- 1
  l <- rep(0, 2 * maxage)
  temp <- rep(1, maxage) - colSums(Umat)
  M <- diag(temp)
  N <- solve(diag(maxage) - Umat)
  eta1 <- rep(1, maxage) %*% N
  B <- M %*% N
  etadagger <- eta1 %*% B
  eta0 <- eta1 %*% e0
  H_new <- etadagger[1] / eta0

  return(as.vector(H_new))
}
