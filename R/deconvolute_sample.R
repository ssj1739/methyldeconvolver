#' deconvolution step
#'
#' @param sample.pat.path 
#' @param reference 
#' @param verbose logical, should messages be outputted to console? Default is False.
#' @param num_of_inits numeric - how many random prior initializations to set for the EM approach. Default is 1 (uniform prior).
#' @param max_iter numeric - how many iterations of EM should be the maximum (Default is 100 - convergence usually achieved before this)
#' @param n_threads numeric - how many cores can be used to parallelize? See details.
#' 
#' @details For parallelization, the \link[pbapply] package is used to speed up computations.
#' Each initialization is sent to a different core for processing, but the EM itself occurs on a single core.
#' To speed up processing with a single core, fewer initialization can be set.
#'
#' @return list of alpha from each initialization.
#' @export
#'
#' @examples
#' \dontrun{
#' ref = learn_reference(marker.file = "my_marker_file.txt", pat.dir = "directory/to/patfiles/")
#' deconvolute_sample_weighted(sample.pat.path = "sample_to_deconvolute.pat.gz", reference = ref)
#' }
deconvolute_sample <- function(sample.pat.path, 
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
    on.exit(pbapply::closepb(pb))
  
  if(verbose) message("Performing EM")
  
  ### Set different random uniform priors based on the number of initializations specified by user 
  num_celltypes <- length(unique(reference$marker$target))
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
  updated_alphas <- pbapply::pblapply(alpha.inits, function(alpha){
    i.iter = 1
    if(max_iter <= i.iter){
      stop("Number of maximum iterations must be greater than 1")
    }
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
      cs = colSums(phi * sample.pat$nobs[omp$overlaps@from])
      new.alpha <- cs / sum(cs)
      
      # Increment iteration
      i.iter = i.iter + 1
      
      # Check our threshold of mad
      mad[i.iter] = mean(abs(alpha.old - new.alpha))/mean(new.alpha)
      alpha = new.alpha
      all.alphas[[i.iter]] <- alpha
    }
    return(list(last_alpha = alpha, iter_mad = mad))
  }, cl = n_threads)
  
  return(updated_alphas)
}
