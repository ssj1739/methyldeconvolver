#' deconvolution step
#'
#' @param sample_pat Path to sample PAT file (output from wgbstools)
#' @param reference Output from [learn_reference()] function.
#' @param filter_reads_with_n_cpgs integer - pass to filter.length, filter reads from PAT file with less than n CpG sites. Default 3.
#' @param quiet logical - should output be delivered silently? Default is FALSE.
#' @param num_of_inits numeric - how many random prior initializations to set 
#' for the EM approach. Default is 1 (uniform prior).
#' @param retain_alphas logical - should alphas from each iteration be saved? Default is FALSE.
#' @param output_format Default is "all". Change to "simple" for only named vector 
#' of proportion estimates
#' @param max_iter numeric - how many iterations of EM should be the maximum 
#' (Default is 100 - convergence usually achieved before this)
#' @param n_threads numeric - how many threads/cores can be used to parallelize? 
#' See details.
#' @param calculate_confidence_int numeric - level at which confidence interval should be calculated (e.g. 0.95). Defaults to NA, which skips the procedure entirely.
#' 
#' @details The [pbapply][pbapply::pbapply] function is used to speed up 
#' computations.
#' Each initialization is sent to a different core for processing, but the EM 
#' itself occurs on a single core.
#' To speed up processing with a single core, 
#' fewer initialization can be set.
#'
#' @return list of alpha from each initialization.
#' @export
#'
#' @examples
#' \dontrun{
#' ref = learn_reference(marker.file = "my_marker_file.txt", pat.dir = "directory/to/patfiles/")
#' deconvolute_sample_weighted(sample_pat = "sample_to_deconvolute.pat.gz", reference = ref)
#' }
deconvolute_sample <- function(sample_pat, 
                               reference,
                               filter_reads_with_n_cpgs = 3,
                               quiet = F, 
                               retain_alphas = F,
                               output_format = "simple",
                               num_of_inits = 10, 
                               max_iter = 100,
                               use.empirical = F,
                               calculate_confidence_int = NA,
                               n_threads = 1){

  # Load required packages
  require(dplyr)
  require(parallel)
  require(stats)
  
  
  # TODO: VALIDATE ALL PARAMS
  if(max_iter <= 1 & !is.integer(max_iter))
    stop("\nNumber of maximum iterations must be integer greater than 1")
  if(num_of_inits < 1 & !is.integer(num_of_inits))
    stop("\nNumber of initializations must be an integer greater than 0")
  if(!is.data.frame(sample_pat)){
    if(!file.exists(sample_pat))
      stop("\nNo such PAT file or dataframe exists.")
  }
  
  # Check available cores (n), do not exceed n-1 cores.
  if(n_threads > parallel::detectCores()){
    if(!quiet) message(paste0("Max number of threads detected is ", 
                              parallel::detectCores(), 
                              ". Using one less than that number instead."))
    n_threads <- parallel::detectCores()-1
  }
  
  # Default to n_threads=1 if platform is not unix
  if(.Platform$OS.type!="unix"){
    n_threads <- 1
    if(!quiet) message("To use multi-core capability of pbapply,
                       UNIX-like systems like Linux or macOS must be used. Changing number of threads to 1.")
  }
  
  if(require(pbapply)){
    xapply <- function(...) pbapply::pblapply(..., cl = n_threads)
  }else{
    if(!quiet) warning("pbapply not installed or available. Using base lapply instead.")
    xapply <- function(...) lapply(...)
  }
  
  if(!quiet) message("Reading and aligning PAT file to marker reference.")
  # Read and validate PAT file
  if(!is.data.frame(sample_pat)){
    sample_pat <- try({
        read_pat(path = sample_pat, filter.noninf = T, filter.length = filter_reads_with_n_cpgs, verbose = !quiet)
    }, silent = T)
    if(!is.data.frame(sample_pat)){
      stop("Unable to read sample pat file. Please provide valid data.frame or path to pat file.")
    }
  }
  # Align markers to reads on PAT file
  omp <- overlap_marker_pat(pat = sample_pat, marker = reference$marker, n_threads = n_threads, split_reads = F)
  
  # Calculate coverage at each marker region
  omp$marker.gr$num_of_reads_per_marker <- sapply(1:length(omp$marker.gr), function(i){
    return(sum(omp$pat.gr$nobs[omp$overlaps@from[omp$overlaps@to==i]]))
  })
  
  if(!quiet) message("Computing psi matrix...")
  
  # Compute psi value for each read, and create new psi matrix
  #psi.mat <- matrix(nrow = length(omp$overlaps), ncol = length(reference$beta_celltype_fits), dimnames = list(c(), names(reference$beta_celltype_fits)))
 
  # Add progress bar for interactive runs
  # TODO: Some regions in the PAT file have multiple overlaps with marker regions, for one of two reasons:
  # - Long reads in the PAT file can span multiple marker regions
  # - The same exact genomic coordinates can be informative markers for multiple cell types
  # A psi value should be computed using   

  res <- xapply(seq_along(omp$overlaps), function(i){
    # Calculate r
    r.vec = encode_binary(omp$pat.gr$read[omp$overlaps@from[i]])
    
    # Extract shape1, shape2, and beta.f for each marker that's associated with the ith read
    # Length of each of these vectors is the # of cell types
    if(use.empirical){
      shape1 = sapply(reference$beta_celltype_fits, function(x) return(x$shape1.emp[omp$overlaps@to[i]]))
      shape2 = sapply(reference$beta_celltype_fits, function(x) return(x$shape2.emp[omp$overlaps@to[i]]))
      beta.f = sapply(reference$beta_celltype_fits, function(x) return(x$beta.f.emp[omp$overlaps@to[i]]))
    }else{
      shape1 = sapply(reference$beta_celltype_fits, function(x) return(x$shape1[omp$overlaps@to[i]]))
      shape2 = sapply(reference$beta_celltype_fits, function(x) return(x$shape2[omp$overlaps@to[i]]))
      beta.f = sapply(reference$beta_celltype_fits, function(x) return(x$beta.f[omp$overlaps@to[i]]))
    }
    
    # TODO: TEMPORARY FIX - replace beta.f values of 0 with small number
    beta.f[beta.f == 0] <- 1e-100
    
    psi.vec = unlist(sapply(reference$beta_celltype_fits, function(x) return(x$psi.init[omp$overlaps@to[i]])))
    psi.vec <- ifelse(is.na(beta.f), 0, 1)
    meth.frac = sum(r.vec) / length(r.vec)
    
    # For all cell types, iterate through each cell type c
    for(c in seq_along(psi.vec)){
      if(psi.vec[c]!=0){
        #P = dbinom(x = sum(r.vec), size = length(r.vec), prob = meth.frac[[c]])

        # For each CpG site r.i in read r (represented in r.vec)
        for(r.i in r.vec){
          # Compute P from the beta function
          P = beta(r.i + shape1[[c]], 1-r.i+shape2[[c]])/beta.f[[c]]
          # Update the psi value for each cell type by multiplying by P
          psi.vec[c] <- psi.vec[c] * P
        }
        psi.vec[c] <- psi.vec[c] * P
      }
    }
    # Output updated psi value into psi.mat
  #  psi.mat[i,] <- psi.vec
    return(list(psi.vec = psi.vec, shape1 = shape1, shape2 = shape2, beta.f = beta.f, meth.frac = meth.frac))
  })
  
  psi.mat <- lapply(res, '[[', "psi.vec") %>% 
    dplyr::bind_rows() %>% 
    as.matrix()
  
  # if(!quiet)
  #   pb <- pbapply::startpb(min = 0, max = length(omp$overlaps))
  # for(i in seq_along(omp$overlaps)){
  #   # Calculate r
  #   r.vec = encode_binary(omp$pat.gr$read[omp$overlaps@from[i]])
  #   
  #   # Extract shape1, shape2, and beta.f for each marker that's associated with the ith read
  #   # Length of each of these vectors is the # of cell types
  # 
  #   # shape1 = sapply(reference$beta_celltype_fits, function(x) return(x$shape1[x$marker.index==omp$overlaps@to[i]]))
  #   # shape2 = sapply(reference$beta_celltype_fits, function(x) return(x$shape2[x$marker.index==omp$overlaps@to[i]]))
  #   # beta.f = sapply(reference$beta_celltype_fits, function(x) return(x$beta.f[x$marker.index==omp$overlaps@to[i]]))
  #   # psi.vec = unlist(sapply(reference$beta_celltype_fits, function(x) return(x$psi.init[x$marker.index==omp$overlaps@to[i]])))
  #   shape1 = sapply(reference$beta_celltype_fits, function(x) return(x$shape1[omp$overlaps@to[i]]))
  #   shape2 = sapply(reference$beta_celltype_fits, function(x) return(x$shape2[omp$overlaps@to[i]]))
  #   beta.f = sapply(reference$beta_celltype_fits, function(x) return(x$beta.f[omp$overlaps@to[i]]))
  #   psi.vec = unlist(sapply(reference$beta_celltype_fits, function(x) return(x$psi.init[omp$overlaps@to[i]])))
  #   meth.frac = sum(r.vec) / length(r.vec)
  #   
  #   # For all cell types, iterate through each cell type c
  #   for(c in seq_along(psi.vec)){
  #     if(psi.vec[c]!=0){
  #       #P = dbinom(x = sum(r.vec), size = length(r.vec), prob = meth.frac[[c]])
  #       
  #       # For each CpG site r.i in read r (represented in r.vec)
  #       for(r.i in r.vec){
  #         # Compute P from the beta function
  #         P = beta(r.i + shape1[[c]], 1-r.i+shape2[[c]])/beta.f[[c]]
  #         # Update the psi value for each cell type by multiplying by P
  #         psi.vec[c] <- psi.vec[c] * P
  #       }
  #       psi.vec[c] <- psi.vec[c] * P
  #     }
  #   }
  #   # Output updated psi value into psi.mat
  #   psi.mat[i,] <- psi.vec
  #   
  #   if(!quiet)
  #     pbapply::setpb(pb = pb, value = i)
  # }
  # if(!quiet)
  #   on.exit(pbapply::closepb(pb))
  
  if(!quiet) message("\nPerforming EM")
  
  ### Set different random uniform priors based on the number of initializations specified by user 
  # Infer number of reference cell types based on number of reference parameter sets learned from cell types
  num_celltypes <- length(reference$beta_celltype_fits)
  # Set first prior to uniform prior
  alpha.inits <- lapply(1:(num_of_inits+1), function(x){
    alpha <- stats::runif(num_celltypes)
    alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) )) # Normalize priors to 1
    names(alpha) <- names(reference$beta_celltype_fits)
    return(alpha)
  })
  alpha.inits[[1]] <- rep(1, num_celltypes) / num_celltypes
  names(alpha.inits[[1]]) <- names(reference$beta_celltype_fits)
  
  
  ### Update alpha to calculate cell type proportions
  updated_alphas <- xapply(alpha.inits, function(alpha){
    i.iter = 1
    mad = vector()
    ll = vector()
    mad[1] <- 1
    all.alphas <- list()
    all.alphas[[i.iter]] <- alpha
    tol = 1e-3
    while(i.iter < max_iter & mad[i.iter] > tol){
      # Save old alpha
      alpha.old = alpha
    
      transpose_prod <- sweep(psi.mat, MARGIN = 2, alpha.old, '*')
      phi <- transpose_prod/rowSums(transpose_prod)

      # Calculate new alpha
      cs = colSums(phi * sample_pat$nobs[omp$overlaps@from])
      new.alpha <- cs / sum(cs)
      
      # Store new alpha
      alpha = new.alpha
      all.alphas[[i.iter]] <- alpha
      
      # Calculate log-likelihood
      ll[i.iter] = likelihood_fun(psi = psi.mat, alpha = new.alpha, epsilon = 1e-99)
      
      # Increment iteration
      i.iter = i.iter + 1
      
      # Check our threshold of mad
      mad[i.iter] = mean(abs(alpha.old - new.alpha))/mean(new.alpha)
    }
    
    output = list(last_alpha = alpha, iter_mad = mad, loglik = ll)
    
    if(retain_alphas)
      output$all_alphas <- all.alphas
    
    return(output)
  })
  
  output <- list(result = updated_alphas)
  
  max_ll <- which.max(sapply(updated_alphas, function(x) max(x[["loglik"]])))
  
  if(!is.na(calculate_confidence_int)){
    if(!is.numeric(calculate_confidence_int) | calculate_confidence_int > 1 | calculate_confidence_int < 0)
      message("Confidence interval specified is invalid. Please specify a numeric confidence interval between 0 and 1.")
    else{
      alpha = calculate_confidence_int
      n.boots = 10 # TODO: Set as parameter
      boot_output <- xapply(1:n.boots, function(x){
        ii = sample(1:nrow(sample_pat), replace = T)
        boot_pat <- sample_pat[ii,]
        
        boot_res <- deconvolute_sample(sample_pat = boot_pat, reference = reference, quiet = T, num_of_inits = 10, n_threads = 1, retain_alphas = F, output_format = 'simple', calculate_confidence_int = NA)
        return(boot_res)
      }, cl = n_threads)
      boot_means <- boot_output %>% dplyr::bind_rows() %>% dplyr::summarize(dplyr::across(dplyr::everything(), mean)) %>% unlist()
      boot_sd <- boot_output %>% dplyr::bind_rows() %>% dplyr::summarize(dplyr::across(dplyr::everything(), stats::sd)) %>% unlist()
      upper_CI <- boot_means + alpha*(boot_sd/sqrt(n.boots))
      lower_CI <- boot_means - alpha*(boot_sd/sqrt(n.boots))
      output$boot_CI <- data.frame(
        boot_means = boot_means,
        boot_sd = boot_sd,
        upper_CI = upper_CI,
        lower_CI = lower_CI
      )
      }
  }
  output$omp <- omp
  if(output_format=="simple")
    return(output$result[[max_ll]]$last_alpha)
  else{
    output$best_result <- output$result[[max_ll]]$last_alpha
    return(output)
  }
}
