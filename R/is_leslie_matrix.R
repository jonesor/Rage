#' Determine if a matrix is a Leslie matrix population model
#'
#' Checks if a given matrix is a Leslie matrix.
#'
#' @details
#' A formal Leslie matrix satisfies the following conditions:
#' * All elements are non-negative.
#' * The subdiagonal elements (survivals), excluding the last column, are all
#'   between 0 and 1.
#' * The sum of the elements in the first row (fecundities) is positive (unless
#'   `includes_mat_F = FALSE`).
#' * The upper triangle, excluding the first row, contains only 0s.
#' * The diagonal, excluding the top-left corner (`A[1,1]`), contains only 0s.
#' * The lower triangle, excluding the subdiagonal, contains only 0s.
#'
#' **Plus-group extension.** In some models the oldest age class is an open-
#' ended *plus group* representing individuals that have survived beyond the
#' maximum explicitly modelled age. These individuals remain in the last class
#' with probability `A[n,n] > 0`. This function accepts such matrices by also
#' exempting the bottom-right corner (`A[n,n]`) from the diagonal-zero
#' requirement. When a plus-group matrix is passed to `leslie_collapse`, the
#' `A[n,n]` term is correctly carried through the internal `leslie_expand` step.
#'
#' @param A Matrix to be tested
#' @param includes_mat_F A logical argument (default `TRUE`) indicating whether
#'   A is expected to include fecundity. The idea here is that A may not include
#'   fertility, but could still be a valid Leslie matrix if fertility was truly
#'   measured to be 0, or if fertility was not measured at all. Thus, this
#'   argument relaxes the test for the first row of A summing to a positive
#'   value.
#'
#' @return A logical value indicating whether the matrix is a Leslie matrix or
#' not
#'
#' @examples
#' # Standard Leslie matrix
#' A <- matrix(c(
#'   0.1, 1.2, 1.1,
#'   0.1, 0.0, 0.0,
#'   0.0, 0.2, 0.0
#' ), nrow = 3, byrow = TRUE)
#' is_leslie_matrix(A) # true
#'
#' # Plus-group Leslie matrix: A[n,n] > 0 is permitted
#' A_plus <- matrix(c(
#'   0.1, 1.2, 1.1,
#'   0.1, 0.0, 0.0,
#'   0.0, 0.2, 0.3
#' ), nrow = 3, byrow = TRUE)
#' is_leslie_matrix(A_plus) # true
#'
#' # Not a Leslie matrix: non-zero off-diagonal elements
#' A <- matrix(c(
#'   0.1, 1.2, 1.1,
#'   0.1, 0.2, 0.1,
#'   0.2, 0.3, 0.3
#' ), nrow = 3, byrow = TRUE)
#' is_leslie_matrix(A) # false
#'
#' # leslie_mpm1 is a plus-group Leslie matrix (A[n,n] > 0)
#' data(leslie_mpm1)
#' A <- leslie_mpm1$matU + leslie_mpm1$matF
#' is_leslie_matrix(A) # true
#'
#' @family transformation
#' @export is_leslie_matrix
#' @author Owen Jones <jones@biology.sdu.dk>
is_leslie_matrix <- function(A, includes_mat_F = TRUE) {
  # Validation of input
  # A must be square matrix and at least dimension 2
  if (!is.matrix(A) || nrow(A) != ncol(A) || nrow(A) < 2) {
    stop("A must be a square matrix with at least dimension 2")
  }

  n <- nrow(A)
  # non-negative elements
  test1 <- !any(A < 0)

  # Subdiagonal - all between 0 and 1
  if (n > 2) {
    sd <- diag(A[-1, -ncol(A)])
  }
  if (n == 2) {
    sd <- A[2, 1]
  }
  test2 <- !any(sd < 0, sd > 1)

  # top row - sum is positive
  if (includes_mat_F) {
    test3 <- sum(A[1, ]) > 0
  } else {
    test3 <- TRUE
  }

  upper_tri <- upper.tri(A)
  upper_tri[1, ] <- FALSE

  # zeros in upper triangle, not including first row
  test4 <- all(A[upper_tri] == 0)

  # zeros in the diagonal, except [1,1] and [n,n] (plus-group corner)
  diagIndicator <- diag(TRUE, n)
  diagIndicator[1, 1] <- FALSE
  diagIndicator[n, n] <- FALSE
  A[diagIndicator]
  test5 <- all(A[diagIndicator] == 0)

  # zeros in lower triangle (not including the survival sub-diagonal)
  lower_tri_1 <- lower.tri(A)
  sd <- data.frame(
    r = 2:nrow(lower_tri_1),
    cc = 1:(nrow(lower_tri_1) - 1)
  )
  for (i in seq_len(nrow(sd))) {
    lower_tri_1[sd$r[i], sd$cc[i]] <- FALSE
  }
  test6 <- all(A[lower_tri_1] == 0)

  return(all(test1, test2, test3, test4, test5, test6))
}
