#' deconvolution step
#'
#' @param sample.pat.path 
#' @param reference 
#' @param num_of_inits numeric - how many random prior initializations to set for the EM approach. Default is 1 (uniform prior).
#'
#' @return list of alpha from each initialization.
#' @export
#'
#' @examples
deconvolute_sample_weighted <- function(sample.pat.path, 
                               reference, 
                               verbose = F, 
                               num_of_inits = 10, max_iter = 100,
                               n_threads = 1){
  # LOAD LIBRARIES
  require(pbapply)
  
  # TODO: VALIDATE PARAMS
  
  if(verbose) message("Reading and aligning PAT file to marker reference.")
  # Read and validate PAT file
  sample.pat <- read_pat(path = sample.pat.path)
  # Align markers to reads on PAT file
  omp <- overlap_marker_pat(pat = sample.pat, marker = reference$marker)
  
  if(verbose) message("Computing psi matrix...")
  
  # Compute psi value for each read, and create new psi matrix
  psi.mat <- matrix(nrow = length(omp$overlaps), ncol = length(reference$beta_celltype_fits), dimnames = list(c(), names(reference$beta_celltype_fits)))
 
  # Compute y value, which will be used in weighting
  y.mat <- matrix(nrow = length(omp$overlaps), ncol = 2, dimnames = list(c(), c("y", "target.ind")))
  #y.mat <- data.frame(y = numeric(length(omp$overlaps)), target = character(length(omp$overlaps)))
  
  if(verbose)
    pb <- pbapply::startpb(min = 0, max = length(omp$overlaps))
  for(i in 1:length(omp$overlaps)){
    # Calculate r
    r.vec = encode_binary(omp$pat.gr$read[omp$overlaps@from[i]])
    
    # Extract shape1, shape2, and beta.f for each marker that's associated with the ith read
    # Length of each of these vectors is the # of cell types
    shape1 = sapply(reference$beta_celltype_fits, function(x) return(x$shape1[x$marker.index==omp$overlaps@to[i]]))
    shape2 = sapply(reference$beta_celltype_fits, function(x) return(x$shape2[x$marker.index==omp$overlaps@to[i]]))
    beta.f = sapply(reference$beta_celltype_fits, function(x) return(x$beta.f[x$marker.index==omp$overlaps@to[i]]))
    psi.vec = unlist(sapply(reference$beta_celltype_fits, function(x) return(x$psi.init[x$marker.index==omp$overlaps@to[i]])))
    meth.frac = sum(r.vec) / length(r.vec)
    
    # Compute weights and populate y.mat
    mu = sapply(reference$beta_celltype_fits, function(x) return(x$mu[x$marker.index==omp$overlaps@to[i]]))
    sigma = sapply(reference$beta_celltype_fits, function(x) return(x$sigma[x$marker.index==omp$overlaps@to[i]]))
    
    n.sigma = 3 ## TODO: PARAMETERIZE THIS
      
    target <- omp$marker.gr$target[omp$overlaps@to[i]]
    upper.bound <- mu[target] + n.sigma * sigma[target]
    lower.bound <- mu[target] - n.sigma * sigma[target]
      
    y.mat[i,1] <- as.numeric(meth.frac < upper.bound & meth.frac > lower.bound)
    y.mat[i,2] <- which(names(mu) == target)
    
    
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
    
    if(verbose)
      pbapply::setpb(pb = pb, value = i)
  }
  if(verbose)
    on.exit(closepb(pb))
  
  if(verbose) message("Performing EM")
  
  ### Set different random uniform priors based on the number of initializations specified by user 
  num_celltypes <- length(unique(reference$marker$target))
  # Set first prior to uniform prior
  alpha.inits <- lapply(1:(num_of_inits+1), function(x){
    alpha <- runif(num_celltypes)
    alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) )) # Normalize priors to 1
    names(alpha) <- names(reference$beta_celltype_fits)
    return(alpha)
  })
  alpha.inits[[1]] <- rep(1, num_celltypes) / num_celltypes
  names(alpha.inits[[1]]) <- names(reference$beta_celltype_fits)
  
  
  ### Update alpha to calculate cell type proportions
  updated_alphas <- pbapply::pblapply(alpha.inits, function(alpha){
    i.iter = 1
    if(max_iter <= i.iter){
      stop("Number of maximum iterations must be greater than 1")
    }
    mad.0 = mad.1 = vector()
    mad.1[1] = mad.0[1] = 1
    
    all.alphas.1 <- list()
    all.alphas.0 <- list()
    all.alphas <- list()
    
    all.alphas.1[[i.iter]] <- alpha
    all.alphas.0[[i.iter]] <- alpha
    all.alphas[[i.iter]] <- alpha
    
    tol = 1e-3
    alpha.1 = alpha
    alpha.0 = alpha
    while(i.iter < max_iter & mad.1[i.iter] > tol & mad.0[i.iter] > tol){
      # Save old alpha
      alpha.old = alpha

      alpha.old.1 = alpha.1
      alpha.old.0 = alpha.0
        
      # Learn weights to apply to each read
      learn.zeta.mat <- data.frame(y = y.mat[,1], alpha = alpha.old[y.mat[,2]])
      zeta.fit <- glm(y ~ alpha, data = learn.zeta.mat, family="binomial")
        
      omega.vec.1 <- predict(zeta.fit, newdata = data.frame(alpha = alpha.old.1), type = "response")
      omega.vec.0 <- predict(zeta.fit, newdata = data.frame(alpha = alpha.old.0), type = "response")
      
      transpose_prod.1 <- sweep(psi.mat[y.mat[,1]==1,], MARGIN = 2, alpha.old.1 * omega.vec.1, '*')
      transpose_prod.0 <- sweep(psi.mat[y.mat[,1]==0,], MARGIN = 2, alpha.old.0 * (1-omega.vec.0), '*')
        
        # Calculate new phi
      #  phi <- transpose_prod/rowSums(transpose_prod)
       phi.1 <- transpose_prod.1/rowSums(transpose_prod.1)
       phi.0 <- transpose_prod.0/rowSums(transpose_prod.0)

       # Calculate new alphas

       cs.1 = colSums(phi.1 * sample.pat$nobs[y.mat[,1]==1, omp$overlaps@from])
       new.alpha.1 = cs.1 / sum(cs.1)
       
       cs.0 = colSums(phi.0 * sample.pat$nobs[y.mat[,1]==0,, omp$overlaps@from])

       new.alpha.0 = cs.0 / sum(cs.0)
     
      omega.vec.1 <- predict(zeta.fit, newdata = data.frame(alpha = new.alpha.1), type = "response")
      omega.vec.0 <- predict(zeta.fit, newdata = data.frame(alpha = new.alpha.0), type = "response")
      
      alpha.avg = ((new.alpha.1 * omega.vec.1) + (new.alpha.0 * (1-omega.vec.0)))/(omega.vec.1 + (1-omega.vec.0))
      new.alpha = alpha.avg/sum(alpha.avg)
        
      i.iter = i.iter + 1
      
      # Check our threshold of mad
      mad.1[i.iter] = mean(abs(alpha.old.1 - new.alpha.1))/mean(new.alpha.1)
      mad.0[i.iter] = mean(abs(alpha.old.0 - new.alpha.0))/mean(new.alpha.0)      
     
      alpha = new.alpha
      
      all.alphas.1[[i.iter]] <- alpha.1
      all.alphas.0[[i.iter]] <- alpha.0
      all.alphas[[i.iter]] <- alpha
    
    }
    return(list(last_alpha = alpha, iter_mad = mad))
  }, cl = n_threads)
  
  return(updated_alphas)
}
