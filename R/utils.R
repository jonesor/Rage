
#' @noRd
checkValidMat <- function(M,
                          fail_all_na = TRUE,
                          fail_any_na = TRUE,
                          warn_all_zero = TRUE,
                          warn_surv_issue = FALSE) {
  
  mn <- deparse(substitute(M)) # name of object passed to M
  
  if (!is.matrix(M) || !is.numeric(M) || (nrow(M) != ncol(M))) {
    stop("Argument ", mn, " must be a square numeric matrix", call. = FALSE)
  }
  if (fail_all_na && all(is.na(M))) {
    stop("Argument ", mn, " contains only missing values (i.e. all <NA>)",
         call. = FALSE)
  }
  if (fail_any_na && any(is.na(M))) {
    stop("Argument ", mn, " contains missing values (i.e. <NA>)", call. = FALSE)
  }
  if (warn_all_zero && all(M == 0)) {
    warning("All elements of ", mn, " are zero", call. = FALSE)
  }
  if (warn_surv_issue && any(colSums(M) > 1)) {
    warning("Argument ", mn, " has at least one stage-specific survival",
            " probability > 1", call. = FALSE)
  }
}

#' @noRd
checkMatchingStageNames <- function(M) {
  if (!identical(rownames(M), colnames(M))) {
    stop("When naming lifestages, both rows and columns must be named ",
         "and their names must be identical")
  }
}

#' @noRd
checkValidStartLife <- function(s, M, start_vec = FALSE) {
  
  if (start_vec) {
    # check that abundances of all stages have been set if passing a population vector
    if ((length(s) > 1 & length(s) != ncol(M)) ||
        (is.numeric(s) && length(s) == 1 && !(s %in% seq_len(nrow(M)))) ||
        (is.character(s) && length(s) == 1 && !(s %in% unique(unlist(dimnames(M)))))) {
      stop("Argument 'start' must be an integer within 1:nrow(matU), ",
           "a character matching a stage class name in dimnames(matU), ",
           "or an integer vector of starting abundances of length ncol(matU)",
           call. = FALSE)
    }
  } else {
    # check that start is a single value, either an index or named life stage
    if ( length(s) != 1 || !(s %in% seq_len(nrow(M))) && !(s %in% unique(unlist(dimnames(M))))) {
      stop("Argument 'start' must be an integer within 1:nrow(matU), ",
           "or a character matching a stage class name in dimnames(matU)",
           call. = FALSE)
    }
  }
}


#' @noRd
colSums2 <- function(mat) {
# like colSums(x, na.rm = TRUE), except that column with only NAs will return NA
#  rather than 0
  apply(mat, 2, function(x) {
    ifelse(all(is.na(x)), NA_real_, sum(x, na.rm = TRUE))
  })
}


### utilities to calculate mean matrix from list of matrices
#' @noRd
meanMat <- function(x, na.rm = FALSE) {
  n_row <- vapply(x, nrow, numeric(1))
  n_col <- vapply(x, ncol, numeric(1))
  if (length(unique(n_row)) != 1 | length(unique(n_col)) != 1) {
    stop("All matrices in list must be of same dimension")
  }
  if (na.rm) x <- lapply(x, zero_NA)
  n <- length(x)
  return(Reduce("+", x) / n)
}


#' @noRd
zero_NA <- function(m) {
  m[is.na(m)] <- 0
  return(m)
}


#' @noRd
area_under_curve <- function(x, y) {
  delta_x <- diff(x)
  if (any(delta_x <= 0)) stop("area_under_curve: x should be ascending")
  rect_heights <- (y[-1] + y[-length(y)]) / 2
  return(sum(delta_x * rect_heights))
}


#' @noRd
#' @importFrom MASS ginv
.matrix_inverse <- function(mat) {
  mat_inv <- try(solve(mat), silent = TRUE)
  if ("try-error" %in% class(mat_inv)) {
    mat_inv <- MASS::ginv(mat)
  }
  return(mat_inv)
}
