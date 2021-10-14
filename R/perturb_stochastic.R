#' Calculate stochastic elasticities from a time-series of matrix population
#' models and corresponding population vectors
#'
#' Calculate stochastic elasticities given a time-series of matrix population
#' models and corresponding population vectors, using the method described in
#' Haridas et al. (2009).
#' 
#' @param X_t A list of matrix population models
#' @param u_t A list of corresponding population vectors
#' 
#' @return A list of three matrices:
#'   \item{E}{matrix of stochastic elasticities}
#'   \item{E_mu}{matrix of stochastic elasticities to mean transition rates}
#'   \item{E_sigma}{matrix of stochastic elasticities to the variance in
#'   transition rates}
#' 
#' @author Patrick Barks <patrick.barks@@gmail.com>
#' 
#' @family perturbation analysis
#' 
#' @references Haridas, C. V., Tuljapurkar, S., & Coulson, T. 2009. Estimating
#'   stochastic elasticities directly from longitudinal data. Ecology Letters,
#'   12, 806-812. <doi:10.1111/j.1461-0248.2009.01330.x>
#' 
#' @examples 
#' # generate list of random MPMs
#' N <- 20     # number of years
#' s <- 3      # matrix dimension
#' X <- list() # matrix population model at time t
#' u <- list() # population vector at time t
#' 
#' for(t in 1:N) {
#'   X[[t]] <- matrix(runif(s^2), nrow = s, ncol = s)
#' }
#' 
#' # derive corresponding series of population vectors
#' u <- pop_vectors(X)
#' 
#' # calculate stochastic elasticities
#' perturb_stochastic(X, u)
#' 
#' @export perturb_stochastic
perturb_stochastic <- function(X_t, u_t) {
  
  if (length(X_t) != length(u_t)) {
    stop("Arguments X_t and u_t must be of same length.\n", call. = FALSE)
  }
  if (length(unique(vapply(X_t, nrow, 1))) != 1) {
    stop("All elements of X_t must be of same dimension.\n", call. = FALSE)
  }
  if (length(unique(vapply(u_t, length, 1))) != 1) {
    stop("All elements of u_t must be of same length.\n", call. = FALSE)
  }
  checks <- lapply(X_t, checkValidMat)
  
  N <- length(X_t)    # number of time periods
  s <- nrow(X_t[[1]]) # matrix dimension
  
  # mean and mean-standardized matrices
  X_mu <- replicate(N, meanMat(X_t), simplify = FALSE)
  X_std <- lapply(X_t, function(X) X - X_mu[[1]])
  
  # calculate observed lambda at each time t
  lambda_t <- mapply(function(X, u) sum(X %*% u), X_t, u_t)
  
  E <- E_mu <- E_sigma <- matrix(0, s, s)
  I <- diag(s)
  e <- rep(1, s)
  
  # for each matrix element
  for(i in 1:s) {
    for(j in 1:s) {
      
      ## build perterbation mats for given matrix element [i, j]
      C <- lapply(X_t, build_pert_mats, i = i, j = j)
      C_mu <- lapply(X_mu, build_pert_mats, i = i, j = j)
      C_sigma <- lapply(X_std, build_pert_mats, i = i, j = j)
      
      ## calculate e_R
      e_R <- mapply(e_fn, mat = C, vec = u_t, lam = lambda_t)
      e_R_mu <- mapply(e_fn, mat = C_mu, vec = u_t, lam = lambda_t)
      e_R_sigma <- mapply(e_fn, mat = C_sigma, vec = u_t, lam = lambda_t)
      
      ## calculate e_U
      w_t <- w_t_mu <- w_t_sigma <- list()
      w_t[[1]] <- w_t_mu[[1]] <- w_t_sigma[[1]] <- rep(0, s)
      
      # recursion to obtain w_t
      for(t in 1:(N-1)) {
        w_t[[t+1]] <- 1 / lambda_t[t] * (I - u_t[[t+1]] %*% t(e)) %*%
          (C[[t]] %*% u_t[[t]] + X_t[[t]] %*% w_t[[t]])
        
        w_t_mu[[t+1]] <- 1 / lambda_t[t] * (I - u_t[[t+1]] %*% t(e)) %*%
          (C_mu[[t]] %*% u_t[[t]] + X_t[[t]] %*% w_t_mu[[t]])
        
        w_t_sigma[[t+1]] <- 1 / lambda_t[t] * (I - u_t[[t+1]] %*% t(e)) %*%
          (C_sigma[[t]] %*% u_t[[t]] + X_t[[t]] %*% w_t_sigma[[t]])
      }
      
      e_U <- mapply(e_fn, mat = X_t, vec = w_t, lam = lambda_t)
      e_U_mu <- mapply(e_fn, mat = X_t, vec = w_t_mu, lam = lambda_t)
      e_U_sigma <- mapply(e_fn, mat = X_t, vec = w_t_sigma, lam = lambda_t)
      
      ## calculate stochastic elasticities
      E[i, j] <- mean(e_R + e_U)
      E_mu[i, j] <- mean(e_R_mu + e_U_mu)
      E_sigma[i, j] <- mean(e_R_sigma + e_U_sigma)
    }
  }
  
  return(list(E = E, E_mu = E_mu, E_sigma = E_sigma))
}


# function to create perturbation matrices (C_ijt)
build_pert_mats <- function(X, i, j) {
  C <- matrix(0, nrow(X), ncol(X))
  C[i, j] <- X[i, j]
  return(C)
}


# function to calculate e_R and e_U
e_fn <- function(mat, vec, lam) { sum(mat %*% vec) / lam }
