#' Calculate likelihood
#'
#' @param psi matrix of psi values
#' @param alpha cell type proportions (same length as ncol(alpha))
#' @param epsilon pseudo-probability of psi in case psi = 0
#'
#' @return log-likelihood of input alpha value
#' @export
#'
likelihood_fun <- function(psi, alpha, epsilon = 1e-99){
  stopifnot(length(alpha) == ncol(psi))
  
  log_alpha <- log(alpha)
  psi[psi == 0] <- epsilon
  log_psi <- log(psi)
  
  
  ll <- sum(rowSums(log_psi)) + sum(log_alpha)
  
  # sum_mat <- t(t(log_psi) + log_alpha)
  # ll <- sum(rowSums(sum_mat, na.rm = T))
  
  return(ll)
}
