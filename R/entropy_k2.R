#' Calculate Keyfitz entropy for a matrix population model
#'
#' Computes Keyfitz entropy from the U submatrix of a matrix population model.
#'
#' @param Umat A square numeric matrix representing the U submatrix of a matrix
#'   population model. For age-based matrices, the dimension corresponds to the
#'   number of age classes
#' @param type An argument specifying the type of matrix model used in the
#'   calculations. This is necessary because the calculations for age vs.
#'   stage-based matrices are different. Possible options are `age` and `stage`.
#'   The latter is not yet implemented. Defaults to `age`.
#' @param init_distrib For stage-based matrices, the initial cohort distribution across
#'   stages. This should sum to 1. If it does not sum to 1, the function
#'   rescales it to 1. Defaults to the stable stage distribution if no
#'   information is provided.
#' @param max_age For stage-based matrices, the upper age, in units of the projection
#'   interval. Defaults to 1000 if no information is provided.
#' @param n_is_maxage If TRUE, survival p_n is set to zero. Defaults to FALSE.
#'
#' @return Returns a single numeric value representing the Keyfitz entropy
#' for the given matrix. This value quantifies the dispersion of age at death.
#'
#' @author Lotte de Vries <c.devries@@uva.nl>
#' @author Stefano Giaimo <giaimo@@evolbio.mpg.de>
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
#' data(mpm1)
#' entropy_k2(mpm1$matU, type = "stage")
#' @importFrom utils head
#' @export

entropy_k2 <- function(Umat, type = "age", init_distrib = NULL, max_age = NULL, n_is_maxage = FALSE) {
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
  if (type == "age") {
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

  # Function from Stefano Giaimo
  if (type == "stage") {
    matDim <- ncol(Umat)
    if (is.null(init_distrib)) {
      init_distrib <- rep(1 / matDim, matDim)
    }
    if (length(init_distrib) != nrow(Umat)) {
      stop("init_distrib must have a length of nrow(Umat)")
    }
    if (is.null(max_age)) {
      max_age <- 1000
    }
    
    #Rescale the initial stage distribution to sum to 1.
    init_distrib <- init_distrib/sum(init_distrib)
    
    Idmat <- diag(ncol(Umat))
    Nmat <- solve(Idmat - Umat)
    powUmat <- Reduce("%*%",
      replicate(max_age, Umat, simplify = FALSE),
      init = Idmat,
      accumulate = TRUE
    )
    lastpowUmat <- powUmat[[max_age + 1]]
    fac1 <- sapply(powUmat, function(x) sum(x %*% init_distrib))
    fac1 <- fac1[-1] / head(fac1, -1)
    if (n_is_maxage == TRUE) {
      fac1[max_age] <- 0
    }
    fac1 <- 1 - fac1
    fac2 <- sapply(head(powUmat, -1), function(x) {
      sum((Nmat %*% (x - lastpowUmat)) %*% init_distrib)
    })
    entr <- sum(fac1 * fac2) / sum((Nmat %*% (Idmat - lastpowUmat)) %*% init_distrib)
    return(entr)
  }
}
