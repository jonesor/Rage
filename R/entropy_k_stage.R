#' Calculate Keyfitz entropy for a stage-based matrix population model
#'
#' Computes Keyfitz entropy from the U submatrix of a stage-based (Lefkovitch)
#' matrix population model.
#'
#' @param Umat A square numeric matrix representing the U submatrix of a
#'   stage-based (Lefkovitch) matrix population model.
#' @param init_distrib The initial cohort distribution across
#'   stages. This should sum to 1. If it does not sum to 1, the function
#'   rescales it to 1. Defaults to an equal distribution across stages.
#' @param max_age The upper age, in units of the projection interval. Defaults
#'   to 1000 if no information is provided.
#' @param n_is_maxage If TRUE, survival p_n is set to zero. Defaults to FALSE.
#'
#' @return Returns a single numeric value representing the Keyfitz entropy
#' for the given matrix. This value quantifies the dispersion of age at death.
#'
#' @author Stefano Giaimo <giaimo@@evolbio.mpg.de>
#' @author Owen Jones <jones@@biology.sdu.dk>
#'
#' @references Keyfitz, N. 1977. Applied Mathematical Demography. New York:
#'   Wiley.
#'
#'
#' @family life history traits
#' @examples
#' data(mpm1)
#'
#' entropy_k_stage(mpm1$matU)
#' @importFrom utils head
#' @export

entropy_k_stage <- function(Umat, init_distrib = NULL, max_age = NULL, n_is_maxage = FALSE) {
  if (!inherits(Umat, "matrix")) {
    stop("m must be a matrix")
  }

  if (!is.numeric(Umat)) {
    stop("m must be a numeric matrix")
  }
  if (nrow(Umat) != ncol(Umat)) {
    stop("m must be a square matrix")
  }

  # Function from Stefano Giaimo
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

  # Rescale the initial stage distribution to sum to 1.
  init_distrib <- init_distrib / sum(init_distrib)

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
