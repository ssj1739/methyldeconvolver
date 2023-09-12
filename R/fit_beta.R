#' fit_beta
#' Fits a beta distribution to given methylation signals from sample pat file
#' Uses optim method to find mle and mme using empirical mean and variance
#'
#' @param overlaps.list output from overlap_marker_pat
#' @param pseudo add pseudocount to avoid Inf/-Inf log values (default = 1e-7)
#' @param verbose logical - verbose output? Default is FALSE.
#'
#' @return beta params for each marker region using stats::optim
#'
#' @importFrom stats dbeta
#' @importFrom stats var
#' @importFrom stats optim
#' @importFrom stringr str_count
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' fit_beta(overlaps.list)
#' }
fit_beta <- function(overlaps.list, pseudo = 1e-7, verbose = F){
  require(stats)
  require(dplyr)
  require(stringr)
  
  marker.pat.overlaps <- overlaps.list$overlaps
  pat.ranges <- overlaps.list$pat.gr
  
  if(verbose)
    message("Learning beta distribution parameters.")
  
  #loop through each cell/marker pair
  beta_fits <- lapply(1:length(overlaps.list$marker.gr), function(x){
    pat.index <- marker.pat.overlaps@from[marker.pat.overlaps@to == x]
    if(length(pat.index)==0){
      return(data.frame(
        shape1 = NA,
        shape2 = NA,
        beta.f = NA,
        shape1.emp = NA,
        shape2.emp = NA,
        beta.f.emp = NA,
        marker.index = x,
        mu = NA,
        sigma = NA
      ))
    }
    
    pat.subset <- pat.ranges[pat.index]
    nTs = stringr::str_count(pat.subset$read, "T")
    nCs = stringr::str_count(pat.subset$read, "C")
    meth.fraction <-  nCs / (nCs + nTs) #caluclate read level methylation fractions
    rep.meth.fraction <- rep(meth.fraction, pat.subset$nobs) #repeat methylation fraction number of times unique PAT read was observed
    
    # Prevent taking the log of 0 by adding/subtracting a psuedo count
    rep.meth.fraction[rep.meth.fraction==0] <- rep.meth.fraction[rep.meth.fraction==0]+pseudo
    rep.meth.fraction[rep.meth.fraction==1] <- rep.meth.fraction[rep.meth.fraction==1]-pseudo
    
    fit.mle <- fit.emp <- rep(NA, 2)
    beta.f <- beta.f.emp <- NA
    mu = NA
    sigma = NA
    
    # If there are not at least 3 distinct reads that map to a given site,
    # and not at least 2 unique values (not all 0s or all 1s), then consider modeling the marker:
    ### NOTE: Come back to the filtering step here.
    
    if(length(rep.meth.fraction) >= 10){
      #method of moments beta distr
      estBetaParams <- function(mu, var) {
        alpha <- (((1 - mu) / var) - (1 / mu)) * mu ^ 2
        beta <- alpha * (1 / mu - 1)
        if(alpha < 0 | beta < 0){
          
        }
        if(is.infinite(abs(alpha)) | is.nan(alpha)){
          if(mu > 0.5){
            alpha <- 20
            beta <- 0.05
          }else{
            alpha <- 0.05
            beta <- 20
          }
        }
        return(params = c(alpha, beta))
      }
    
      #beta likeihood function to maximize
      beta_likelihood <- function(theta, x){
        N <- length(x)
        alpha.val <- theta[1]
        beta.val <- theta[2]
        beta_ll <- ((alpha.val-1)*sum((log(x)))) + ((beta.val-1)*sum(log(1-x))) - N*log(beta(alpha.val, beta.val))
        return(-beta_ll)
      }
      
      # Empirical mean and variance of for mme fit
      mu <- mean(rep.meth.fraction)
      sigma <- var(rep.meth.fraction)
      
      # Set lower bound on sigma
      if(sigma < 1e-4){
        sigma <- 1e-4
      }
      
      # Fit beta distr by mle using beta likelihood function and optim function with L-BFGS-B method
      fit.mle <- tryCatch({
        res = stats::optim(par = c(0.01, 0.01), 
                         fn = beta_likelihood, 
                         x = rep.meth.fraction) #, 
                         #method = "L-BFGS-B", 
                         #lower = 0.0001, upper = 10000)
        res$par
      }, error = function(e) return(c(NA, NA)))
      # Set heuristics based on calculated mu if unable to learn parameters
      if(all(is.na(fit.mle))){
        if(mu > 0.5){
          fit.mle[1] <- 20
          fit.mle[2] <- 0.05
        }else{
          fit.mle[1] <- 0.05
          fit.mle[2] <- 20
        }
      }
      if(mu > 0.9999){
        fit.mle[1] <- 20
        fit.mle[2] <- 0.05
      }
      if(mu < 0.0001){
        fit.mle[1] <- 0.05
        fit.mle[2] <- 20
      }
      
      # Empirically calculate beta distr using empirical mean and variance of read level methylation fractions
      fit.emp = estBetaParams(mu, sigma)
      
      # Derive beta function value from fitted beta parameters, used in later computations
      #Currently only returning mle version
      beta.f <- beta(fit.mle[1], fit.mle[2])
      beta.f.emp <- beta(fit.emp[1], fit.emp[2])
    }
    #return both mle and mme estimates with meta data
    return(data.frame(shape1 = fit.mle[1], shape2 = fit.mle[2], beta.f = beta.f, 
                      shape1.emp = fit.emp[1], shape2.emp = fit.emp[2], beta.f.emp = beta.f.emp,
                      marker.index = x, mu = mu, sigma = sigma, n = length(rep.meth.fraction)))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-dplyr::starts_with("X")) %>%
    dplyr::mutate(psi.init = as.numeric(!is.na(beta.f)))
  return(beta_fits)
}



