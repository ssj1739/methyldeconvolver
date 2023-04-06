#' fit_beta_old
#'
#'uses fitdistrplus mle and mme methods
#'
#' @param overlaps.list 
#'
#' @return DEPRECATED - list of beta params using fitdistr
#'
#' @examples
#' \dontrun{
#' fit_beta_old(overlaps.list, pf = 0.05)
#' }
fit_beta_old <- function(overlaps.list, pf = 0.05){
  require(fitdistrplus)
  
  marker.pat.overlaps <- overlaps.list$overlaps
  pat.ranges <- overlaps.list$pat.gr
  
  beta_fits <- lapply(1:length(overlaps.list$marker.gr), function(x){
    pat.index <- marker.pat.overlaps@from[marker.pat.overlaps@to == x]
    if(length(pat.index)==0){
      return(data.frame(
        shape1 = NA,
        shape2 = NA,
        beta.f = NA,
        marker.index = x,
        avg.meth = NA
      ))
    }
    # Potentially add filter step to remove <=2 length of pat.index
    
    pat.subset <- pat.ranges[pat.index]
    nTs = stringr::str_count(pat.subset$read, "T")
    nCs = stringr::str_count(pat.subset$read, "C")
    meth.fraction <-  nCs / (nCs + nTs) #calc read level methylation fraction
    rep.meth.fraction <- rep(meth.fraction, pat.subset$nobs) #repeat methylation fraction by number of times PAT read was observed
    
    # If there are not at least 3 distinct reads that map to a given site,
    # and not at least 2 unique values (not all 0s or all 1s), then consider modeling the marker:
    
    #Add/subtract pseudo count
    avg.meth <- mean(rep.meth.fraction, na.rm = T)
    if(avg.meth==0){
      avg.meth = avg.meth + pf
    }
    if(avg.meth==1){
      avg.meth = avg.meth - pf
    }
    fit <- rep(NA, 2) #beta parameters
    beta.f <- NA #beta function calculation for later computation
    if(length(rep.meth.fraction) >= 3 & length(unique(rep.meth.fraction)) >= 2){
      # First try mle:
      fit = try({
        fitdistrplus::fitdist(data = rep.meth.fraction, distr = "beta", method = "mle")
        }, silent = T)
      
      # If mle doesn't work, then try mme:
      if(class(fit)=="try-error"){ 
        fit = fitdistrplus::fitdist(data = rep.meth.fraction, distr = "beta", method = "mme")
      }
      fit <- fit$estimate
      beta.f <- beta(fit[1], fit[2])
    }
    return(data.frame(t(fit), beta.f = beta.f, marker.index = x, avg.meth = avg.meth))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-dplyr::starts_with("X")) %>%
    dplyr::mutate(psi.init = as.numeric(!is.na(avg.meth)))
  return(beta_fits)
}

#' fit_beta
#' Fits a beta distribution to given methylation signals from sample pat file using optim function to maximize beta distr log-likelihood
#'
#' @param overlaps.list output from overlap_marker_pat
#' @param pseudo add pseudocount to avoid Inf/-Inf log values (default = 1e-7)
#'
#' @return beta params for each marker region using stats::optim
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
  
  beta_fits <- lapply(1:length(overlaps.list$marker.gr), function(x){
    pat.index <- marker.pat.overlaps@from[marker.pat.overlaps@to == x]
    if(length(pat.index)==0){
      return(data.frame(
        shape1 = NA,
        shape2 = NA,
        beta.f = NA,
        marker.index = x,
        mu = NA,
        sigma = NA
      ))
    }
    
    pat.subset <- pat.ranges[pat.index]
    nTs = stringr::str_count(pat.subset$read, "T")
    nCs = stringr::str_count(pat.subset$read, "C")
    meth.fraction <-  nCs / (nCs + nTs) #calc read level methylation fraction
    rep.meth.fraction <- rep(meth.fraction, pat.subset$nobs) #repeat methylation fraction by number of times PAT read was observed
    
    # Prevent taking the log of 0 by adding/subtracting pseudo count
    rep.meth.fraction[rep.meth.fraction==0] <- rep.meth.fraction[rep.meth.fraction==0]+pseudo
    rep.meth.fraction[rep.meth.fraction==1] <- rep.meth.fraction[rep.meth.fraction==1]-pseudo
    
    fit <- rep(NA, 2) #beta distr parameters
    beta.f <- NA #beta function for later computation
    mu = NA #mean of beta distribution
    sigma = NA #variance of beta distribution
    
    # If there are not at least 3 distinct PAT reads that map to a given site,
    # and not at least 2 unique values (not all 0s or all 1s), then consider modeling the marker:
    #
    ### NOTE: Come back to the filtering step here or consider specifying as user input parameters
    #
    #Calculate beta distribution log-likelihood using vector of read levelmethylation fractions at each marker region
    if(length(rep.meth.fraction) >= 3){
      beta_likelihood <- function(theta, x){ #define beta log-likelihood
        N <- length(x)
        alpha <- theta[1]
        beta <- theta[2]
        beta_ll <- ((alpha-1)*sum((log(x)))) + ((beta-1)*sum(log(1-x))) - N*log(beta(alpha, beta))
        return(-beta_ll)
      }
      
      #Use optim function to find mle
      res = stats::optim(par = c(0.01, 0.01), fn = beta_likelihood, x = rep.meth.fraction, method = "L-BFGS-B", lower = 0.001, upper = 10000)
      
      fit <- res$par
      beta.f <- beta(fit[1], fit[2])
      mu <- fit[1]/(fit[1]+fit[2])
      sigma <- sqrt((fit[1]*fit[2])/(((fit[1]+fit[2])^2)*(fit[1]+fit[2]+1)))
      
      
      
    }
    #output alpha, beta, beta function, marker index, mean, variance
    return(data.frame(shape1 = fit[1], shape2 = fit[2], beta.f = beta.f, marker.index = x, mu = mu, sigma = sigma))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-dplyr::starts_with("X")) %>%
    dplyr::mutate(psi.init = as.numeric(!is.na(beta.f)))
  return(beta_fits)
}

#' fit_beta_new
#' Fits a beta distribution to given methylation signals from sample pat file
#' Uses optim method to find mle and mme using empirical mean and variance
#'
#' @param overlaps.list output from overlap_marker_pat
#' @param pseudo add pseudocount to avoid Inf/-Inf log values (default = 1e-7)
#'
#' @return beta params for each marker region using stats::optim
#'
#' @export
#' @examples
#' \dontrun{
#' fit_beta_new(overlaps.list)
#' }
fit_beta_new <- function(overlaps.list, pseudo = 1e-7, verbose = F){
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
    beta.f <- NA
    mu = NA
    sigma = NA
    
    # If there are not at least 3 distinct reads that map to a given site,
    # and not at least 2 unique values (not all 0s or all 1s), then consider modeling the marker:
    ### NOTE: Come back to the filtering step here.
    
    if(length(rep.meth.fraction) >= 3){
      #method of moments beta distr
      estBetaParams <- function(mu, var) {
        alpha <- (((1 - mu) / var) - (1 / mu)) * mu ^ 2
        beta <- alpha * (1 / mu - 1)
        return(params = c(alpha, beta))
      }
      
      #beta likeihood function to maximize
      beta_likelihood <- function(theta, x){
        N <- length(x)
        alpha <- theta[1]
        beta <- theta[2]
        beta_ll <- ((alpha-1)*sum((log(x)))) + ((beta-1)*sum(log(1-x))) - N*log(beta(alpha, beta))
        return(-beta_ll)
      }
      
      # Empirical mean and variance of for mme fit
      mu <- mean(rep.meth.fraction)
      sigma <- var(rep.meth.fraction)
      
      # Fit beta distr by mle using beta likelihood function and optim function with L-BFGS-B method
      res = stats::optim(par = c(0.01, 0.01), fn = beta_likelihood, x = rep.meth.fraction, method = "L-BFGS-B", lower = 0.0001, upper = 10000)
      fit.mle <- res$par
      
      # Empirically calculate beta distr using empirical mean and variance of read level methylation fractions
      fit.emp = estBetaParams(mu, sigma)
      
      # Derive beta function value from fitted beta parameters, used in later computations
      #Currently only returning mle version
      beta.f <- beta(fit.mle[1], fit.mle[2])
    }
    #return both mle and mme estimates with meta data
    return(data.frame(shape1 = fit.mle[1], shape2 = fit.mle[2], beta.f = beta.f, 
                      shape1.emp = fit.emp[1], shape2.emp = fit.emp[2],
                      marker.index = x, mu = mu, sigma = sigma))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-dplyr::starts_with("X")) %>%
    dplyr::mutate(psi.init = as.numeric(!is.na(beta.f)))
  return(beta_fits)
}



