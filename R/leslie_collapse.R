#' Aggregate a Leslie matrix
#'
#' Takes a Leslie matrix and aggregates it to a desired dimension m using least
#' squares weights equal to the stable age distribution. The output includes the
#' aggregated matrix, the weight matrix, the original (or expanded) Leslie
#' matrix raised to the k power, the partitioning matrix, the size of the
#' original (or expanded) Leslie matrix, the size of the aggregated matrix, the
#' quotient of the original size divided by the aggregated size, and the
#' effectiveness of aggregation.
#'
#' @param A a Leslie matrix
#' @param m the dimensionality of the desired aggregated matrix
#'
#' @return a list including the following elements:
#' \item{B}{The aggregated matrix}
#' \item{W}{The weight matrix}
#' \item{Ak}{The original (or expanded) Leslie matrix raised to the k power}
#' \item{S}{The partitioning matrix}
#' \item{n}{The size of the original (or expanded) Leslie matrix}
#' \item{m}{The size of the aggregated matrix}
#' \item{k}{The quotient of the original size divided by the aggregated size}
#' \item{EFF}{The effectiveness of aggregation}
#'
#' @examples
#' data(leslie_mpm1)
#' A <- leslie_mpm1$matU + leslie_mpm1$matF
#' leslie_collapse(A, 4)
#'
#' @author Richard A. Hinrichsen <rich@hinrichsenenvironmental.com>
#'
#' @export leslie_collapse
#' @family transformation
#' @importFrom expm %^%
leslie_collapse <- function(A, m) {
  # data input validation
  if(!m%%1==0 || m<1) {stop("m must be a positive integer")}
  if(nrow(A)!=ncol(A)) {stop("A must be a square matrix")}
  
  # first check whether the matrix is a Leslie matrix
  if (!is_leslie(A)) {stop("A must be a Leslie matrix")}
  n <- dim(A)[1]
  if (n %% m != 0) {
    A <- leslie_expand(A, m)
  } # if
  n <- dim(A)[1]
  k <- n / m
  eigA <- eigen(A)
  iii <- (Re(eigA$values) == max(Re(eigA$values)))
  v <- abs(Re(eigA$vectors[, iii]))
  w <- v # weights are equal to the stable age distribution
  Ak <- A %^% k
  e <- rep(1, k)
  Id <- diag(1, m)
  S <- kronecker(Id, t(e))
  W <- diag(w)
  TT <- W %*% t(S) %*% solve(S %*% W %*% t(S))
  B <- S %*% Ak %*% TT
  EFF <- norm(B %*% S %*% sqrt(W), type = "F")^2 / norm(S %*% Ak %*% sqrt(W),
    type = "F"
  )^2
  return(list(B = B, W = W, Ak = Ak, S = S, n = n, m = m, k = k, EFF = EFF))
}

#' Expand a Leslie matrix
#'
#' Expands a Leslie matrix to a desired dimensionality m by adding subdiagonal
#' elements of 1 and diagonal elements from the original matrix. The output is
#' the expanded matrix.
#'
#' @param A a Leslie matrix
#' @param m the desired dimensionality of the expanded matrix
#'
#' @return the expanded matrix
#'
#' @examples
#' A <- matrix(c(
#'   3.8930214, 8.3389426, 5.651921,
#'   0.5670728, 0.0000000, 0.000000,
#'   0.0000000, 0.2194218, 0.000000
#' ), byrow = TRUE, ncol = 3)
#' leslie_expand(A, 2)
#' @author Richard A. Hinrichsen <rich@hinrichsenenvironmental.com>
#'
#' @noRd
leslie_expand <- function(A, m) {
  n <- dim(A)[1]
  A.expanded <- matrix(0, nrow = n * m, ncol = n * m)
  # fill in first row of A.expanded
  for (ii in 1:n) {
    A.expanded[1, ii * m] <- A[1, ii]
  }
  # next fill in the subdiagonal
  for (ii in 2:(n * m)) {
    A.expanded[ii, ii - 1] <- 1
  }
  for (ii in 1:(n - 1)) {
    A.expanded[m * ii + 1, m * ii] <- A[ii + 1, ii]
  }
  return(A.expanded)
}

#' Check whether a matrix is a Leslie matrix
#'
#' Checks whether a matrix is a Leslie matrix by verifying that the matrix
#' satisfies three conditions: (1) the matrix has non-negative elements, (2) the
#' elements on the subdiagonal are between 0 and 1, and (3) the sum of the birth
#' rates is positive.
#'
#' @param A a matrix to be checked
#'
#' @return a logical value indicating whether the matrix is a Leslie matrix
#'
#' @examples
#' A <- matrix(c(0.1, 1.2, 1.1, 
#'               0.0, 0.2, 0.0, 
#'               0.0, 0.0, 0.3), nrow = 3, byrow = TRUE)
#' is_leslie(A) # true
#' A <- matrix(c(0.1, 1.2, 1.1, 
#'               0.1, 0.2, 0.1, 
#'               0.2, 0.3, 0.3), nrow = 3, byrow = TRUE)
#' is_leslie(A) # false?
#' @author Richard A. Hinrichsen <rich@hinrichsenenvironmental.com>
#'
#' @export is_leslie
is_leslie <- function(A) {
  EPS <- 0.00001
  if (is.na(sum(A))) {
    return(FALSE)
  }
  m <- dim(A)[1]
  if (m == 1) {
    return(TRUE)
  }
  birth.rates <- A[1, ]
  if (m > 2) {
    survival.rates <- diag(A[2:m, 1:(m - 1)])
  }
  if (m == 2) {
    survival.rates <- A[2, 1]
  }

  Anew <- matrix(0, nrow = m, ncol = m)
  if (m > 2) Anew[2:m, 1:(m - 1)] <- diag(survival.rates)
  if (m == 2) {
    Anew[2, 1] <- survival.rates
  }
  Anew[1, ] <- birth.rates

  error <- c(A) - c(Anew)
  error.squared <- t(error) %*% error
  test1 <- (sqrt(error.squared) < EPS)
  iii <- (survival.rates >= 0) & (survival.rates <= 1)
  test2 <- (sum(iii) == (m - 1))
  iii <- (birth.rates >= 0)
  test3 <- (sum(iii) == m)
  return(is.logical(test1 & test2 & test3))
}
