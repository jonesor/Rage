
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
checkValidStartLife <- function(s, M) {
  
  if ( length(s) != 1 || !(s %in% seq_len(nrow(M))) ) {
    stop("Argument startLife must be an integer within 1:nrow(matU)",
         call. = FALSE)
  }
}


#' @noRd
colSums2 <- function(mat) {
  apply(mat, 2, function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    } else {
      return(sum(x, na.rm = TRUE))
    }
  })
}

