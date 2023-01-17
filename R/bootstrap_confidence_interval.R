boostrap_confidence_interval <- function(num_boots = 1000, sample_pat, reference, n_threads){
  pbapply::pblapply(1:num_boots, function(x){
    boot_pat <- sample.pat %>% dplyr::slice_sample(n = nrow(sample.pat))
    boot_omp <- overlap_marker_pat(pat = boot_pat, marker = reference$marker)
    
    psi.mat <- matrix(nrow = length(boot_omp$overlaps), ncol = length(boot_reference$beta_celltype_fits), dimnames = list(c(), names(boot_reference$beta_celltype_fits)))
    if(!quiet)
      pb <- pbapply::startpb(min = 0, max = length(boot_omp$overlaps))
    for(i in 1:length(omp$overlaps)){
      # Calculate r
      r.vec = encode_binary(boot_omp$pat.gr$read[boot_omp$overlaps@from[i]])
      
      # Extract shape1, shape2, and beta.f for each marker that's associated with the ith read
      # Length of each of these vectors is the # of cell types
      shape1 = sapply(reference$beta_celltype_fits, function(x) return(x$shape1[x$marker.index==boot_omp$overlaps@to[i]]))
      shape2 = sapply(reference$beta_celltype_fits, function(x) return(x$shape2[x$marker.index==boot_omp$overlaps@to[i]]))
      beta.f = sapply(reference$beta_celltype_fits, function(x) return(x$beta.f[x$marker.index==boot_omp$overlaps@to[i]]))
      psi.vec = unlist(sapply(reference$beta_celltype_fits, function(x) return(x$psi.init[x$marker.index==boot_omp$overlaps@to[i]])))
      meth.frac = sum(r.vec) / length(r.vec)
      
      # For all cell types, iterate through each cell type c
      for(c in 1:length(psi.vec)){
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
      psi.mat[i,] <- psi.vec
      
      if(!quiet)
        pbapply::setpb(pb = pb, value = i)
    }
    if(!quiet)
      on.exit(pbapply::closepb(pb))
    
    # EM Step
    if(!quiet) message("\nPerforming EM")
    
    ### Set different random uniform priors based on the number of initializations specified by user 
    num_celltypes <- length(unique(reference$marker$target))
    # Set first prior to uniform prior
    num_of_inits = 9
    alpha.inits <- lapply(1:(num_of_inits+1), function(x){
      alpha <- stats::runif(num_celltypes)
      alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) )) # Normalize priors to 1
      names(alpha) <- names(reference$beta_celltype_fits)
      return(alpha)
    })
    alpha.inits[[1]] <- rep(1, num_celltypes) / num_celltypes
    names(alpha.inits[[1]]) <- names(reference$beta_celltype_fits)
  
    ### Update alpha to calculate cell type proportions
    updated_alphas <- pbapply::pblapply(alpha.inits, function(alpha){
      i.iter = 1
      mad = vector()
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
        cs = colSums(phi * boot_pat$nobs[boot_omp$overlaps@from])
        new.alpha <- cs / sum(cs)
        
        # Increment iteration
        i.iter = i.iter + 1
        
        # Check our threshold of mad
        mad[i.iter] = mean(abs(alpha.old - new.alpha))/mean(new.alpha)
        alpha = new.alpha
        all.alphas[[i.iter]] <- alpha
      }
      
      ll = likelihood_fun(psi = psi.mat, alpha = alpha)
      
      output = list(last_alpha = alpha, iter_mad = mad, loglik = ll)
      
      if(retain_alphas)
        output$all_alphas <- all.alphas
      
      return(output)
    }, cl = n_threads) 
  }, cl = n_threads)
}