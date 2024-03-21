#' Calculate Keyfitz entropy for a matrix population model
#'
#' This function computes Keyfitz entropy from a transition matrix (the A
#' matrix). The function is designed to work with discrete, age-classified
#' matrix models.
#'
#' @param A A square numeric projection matrix The dimension of the matrix
#'   corresponds to the number of age classes
#' @param type An argument specifying the type of matrix model used in the
#'   calculations. This is necessary because the calculations for age vs.
#'   stage-based matrices are different. Possible options are `age` and `stage`.
#' 
#' @return Returns a single numeric value representing the Keyfitz entropy
#' for the given matrix. This value quantifies the dispersion of age at death.
#'
#' @author Lotte de Vries <c.devries@@uva.nl>
#' @author Owen R. Jones <jones@@biology.sdu.dk>
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
#' # Define a matrix population model for a hypothetical population
#' A <- matrix(c(
#'   0.7, 0.0, 0.0, 0.5,
#'   0.4, 0.8, 0.0, 0.0,
#'   0.0, 0.2, 0.9, 0.0,
#'   0.0, 0.0, 0.2, 0.1
#' ), nrow = 4, byrow = TRUE)
#'
#' # Calculate Keyfitz entropy
#' entropy_k2(A)
#'
#' @export

entropy_k2 <- function(A, type = "age") {
  if (!inherits(A, "matrix")) {
    stop("A must be a matrix")
  }
  if (!is.numeric(A)) {
    stop("A must be a numeric matrix")
  }
  if (nrow(A) != ncol(A)) {
    stop("A must be a square matrix")
  }
if(type == "age"){
  maxage <- dim(A)[1]
  e0 <- rep(0, maxage)
  e0[1] <- 1
  l <- rep(0, 2 * maxage)
  temp <- rep(1, maxage) - colSums(A)
  M <- diag(temp)
  N <- solve(diag(maxage) - A)
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