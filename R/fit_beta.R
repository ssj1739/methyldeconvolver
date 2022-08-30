#' fit_beta_old
#'
#' @param overlaps.list 
#'
#' @return DEPRECATED - list of beta params using fitdistr
#' @export
#'
#' @examples
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
    meth.fraction <-  nCs / (nCs + nTs)
    rep.meth.fraction <- rep(meth.fraction, pat.subset$nobs)
    
    # If there are not at least 3 distinct reads that map to a given site,
    # and not at least 2 unique values (not all 0s or all 1s), then consider modeling the marker:
    
    avg.meth <- mean(rep.meth.fraction, na.rm = T)
    if(avg.meth==0){
      avg.meth = avg.meth + pf
    }
    if(avg.meth==1){
      avg.meth = avg.meth - pf
    }
    fit <- rep(NA, 2)
    beta.f <- NA
    if(length(rep.meth.fraction) >= 3 & length(unique(rep.meth.fraction)) >= 2){
      # First try mle:
      fit = try(fitdistrplus::fitdist(data = rep.meth.fraction, distr = "beta", method = "mle"), silent = T)
      
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
#'
#' @param overlaps.list output from overlap_marker_pat
#' @param epsilon add pseudocount to avoid Inf/-Inf log values
#'
#' @return beta params for each marker region using stats::optim
#' @export
#'
#' @examples
fit_beta <- function(overlaps.list, epsilon = 1e-7){
  require(stats)
  require(dplyr)
  require(stringr)
  
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
        mu = NA,
        sigma = NA
      ))
    }
    # Potentially add filter step to remove <=2 length of pat.index
    
    pat.subset <- pat.ranges[pat.index]
    nTs = stringr::str_count(pat.subset$read, "T")
    nCs = stringr::str_count(pat.subset$read, "C")
    meth.fraction <-  nCs / (nCs + nTs)
    rep.meth.fraction <- rep(meth.fraction, pat.subset$nobs)
    
    # Prevent taking the log of 0:
    rep.meth.fraction[rep.meth.fraction==0] <- rep.meth.fraction[rep.meth.fraction==0]+epsilon
    rep.meth.fraction[rep.meth.fraction==1] <- rep.meth.fraction[rep.meth.fraction==1]-epsilon
    
    fit <- rep(NA, 2)
    beta.f <- NA
    mu = NA
    sigma = NA
    
    # If there are not at least 3 distinct reads that map to a given site,
    # and not at least 2 unique values (not all 0s or all 1s), then consider modeling the marker:
    ### NOTE: Come back to the filtering step here.
    if(length(rep.meth.fraction) >= 3){
      beta_likelihood <- function(theta, x){
        N <- length(x)
        alpha <- theta[1]
        beta <- theta[2]
        beta_ll <- ((alpha-1)*sum((log(x)))) + ((beta-1)*sum(log(1-x))) - N*log(beta(alpha, beta))
        return(-beta_ll)
      }
      
      res = stats::optim(par = c(0.01, 0.01), fn = beta_likelihood, x = rep.meth.fraction, method = "L-BFGS-B", lower = 0.01, upper = 100)
      
      fit <- res$par
      beta.f <- beta(fit[1], fit[2])
      mu <- fit[1]/(fit[1]+fit[2])
      sigma <- sqrt((fit[1]*fit[2])/(((fit[1]+fit[2])^2)*(fit[1]+fit[2]+1)))
      
    }
    return(data.frame(shape1 = fit[1], shape2 = fit[2], beta.f = beta.f, marker.index = x, mu = mu, sigma = sigma))
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-dplyr::starts_with("X")) %>%
    dplyr::mutate(psi.init = as.numeric(!is.na(beta.f)))
  return(beta_fits)
}
