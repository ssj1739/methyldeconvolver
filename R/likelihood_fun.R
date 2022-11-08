#' Calculate likelihood
#'
#' @param psi matrix of psi values
#' @param alpha cell type proportions (same length as ncol(alpha))
#' @param epsilon pseudo-probability of psi in case psi = 0
#'
#' @return
#' @export
#'
#' @examples
likelihood_fun <- function(psi, alpha, epsilon = 1e-7){
  stopifnot(length(alpha) == ncol(psi))
  
  log_alpha <- log(alpha)
  psi[psi == 0] <- epsilon
  log_psi <- log(psi)
  
  sum_mat <- t(t(log_psi) + log_alpha)
  
  return(sum(rowSums(sum_mat)))
}