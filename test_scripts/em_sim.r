library(methyldeconvolveR)
library(pbapply)

###Simulation psuedo code
#input(num_celltypes, num_markers_per_celltype, num_reads)
num_celltypes = 5
num_markers_per_celltype = 25
num_reads = 10000
num_error_sd = 2
#noise = 0.5

#Generate alpha:
alpha <- runif(num_celltypes)
alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
cell <- seq(1,num_celltypes)
true_alpha <- cbind(cell, alpha)
colnames(true_alpha) <- c('cell_idx', 'alpha')

trial_output <- pblapply(1:100, function(trial){
  noise = 0.01*trial
  #Generate marker target indicator variable column:
  sequence <- c((rep(0, num_celltypes-1)),1)
  permutations <- unique(arrangements::permutations(sequence, layout = "list"))
  target <- rep(unlist(permutations), num_markers_per_celltype)
  
  #matrix(target, nrow = num_markers_per_celltype * num_celltypes, ncol = num_celltypes, byrow = T)
  
  #Simulate reference parameters
  total_markers <- num_celltypes*num_markers_per_celltype
  marker_idx <- rep(1:total_markers, each= num_celltypes)
  cell_idx <- rep(1:num_celltypes, total_markers)
  eta <- vector(mode = "numeric",total_markers*num_celltypes)
  rho <- vector(mode = "numeric", total_markers*num_celltypes)
  mu <- vector(mode = "numeric", total_markers*num_celltypes)
  sigma <- vector(mode = "numeric", total_markers*num_celltypes)
  ref <- cbind(marker_idx, cell_idx, target, eta, rho, mu, sigma)
  colnames(ref) <- c('marker_idx', 'cell_idx', 'target', 'eta', 'rho', 'mu', 'sigma')
  
  for (i in 1:nrow(ref)){
    if (ref[i,'target'] == 1){
      ref[i,'eta'] <- runif(1,0.001,2)
      ref[i, 'rho'] <- runif(1,8,20)
    } else {
      ref[i,'eta'] <- runif(1,8,20)
      ref[i, 'rho'] <- runif(1,0.001,2)
    }
  }
  
  # Calculate mu and sigma of weights
  # ref[i, 'mu'] <- ref[i, 'eta']/(ref[i,'eta']+ref[i,'rho'])
  # ref[i, 'sigma'] <- sqrt((ref[i,'eta']*ref[i,'rho'])/(((ref[i,'eta']+ref[i,'rho'])^2)*(ref[i,'eta']+ref[i,'rho']+1)))
  
  for(i in 1:nrow(ref)){
    mi <- ref[i,"marker_idx"]
    ta_ref <- ref[ref[,"marker_idx"]==mi & ref[,"target"]==1,]
    ref[i, 'mu'] <- ta_ref['eta']/(ta_ref['eta']+ta_ref['rho'])
    ref[i, 'sigma'] <- sqrt((ta_ref['eta']*ta_ref['rho'])/(((ta_ref['eta']+ta_ref['rho'])^2)*(ta_ref['eta']+ta_ref['rho']+1)))
  }
  
  # Transform ref into format used in package
  ref_list <- ref %>%
    as.data.frame() %>%
    dplyr::select(-target, shape1=eta, shape2=rho, marker.index = marker_idx, mu = mu, sigma = sigma) %>%
    dplyr::mutate(beta.f = beta(shape1, shape2), psi.init = 1) %>%
    split(cell_idx)
  
  # #Generate alpha:
  # alpha <- runif(num_celltypes)
  # alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
  # cell <- seq(1,num_celltypes)
  # true_alpha <- cbind(cell, alpha)
  # colnames(true_alpha) <- c('cell_idx', 'alpha')
  # 
  #Generate reads
  #read_idx <- 1
  #colnames(read_mat) <- c("marker.index", "read", "nobs", "target")
  sim_pat <- lapply(1:num_reads, function(x){
    cell_sim <- rmultinom(1,1,true_alpha[,"alpha"]) #simulate cell of origin
    cell_origin <- which(cell_sim == 1) 
    ref_sub_target <- subset(ref, cell_idx == cell_origin & target == 1) #subset reference by cell of origin
    ref_sub_noise <- subset(ref, cell_idx == cell_origin & target == 0) #subset reference by cell of origin
    n <- rbinom(1,1,1-noise)
    if (n == 1){
      ref_sample <- ref_sub_target[sample(nrow(ref_sub_target), 1), ] #sample row from reference subset for alignment
    } else {ref_sample <- ref_sub_noise[sample(nrow(ref_sub_noise), 1), ]}
    
    # ref_sub <- subset(ref, cell_idx == cell_origin) #subset reference by cell of origin
    # ref_sample <- ref_sub[sample(nrow(ref_sub), 1), ] #sample row from reference subset for alignment
    read_length <- round(runif(1,4,12)) #simulate read length
    num_obs <- rpois(1,10) #simulate number of observations of unique read
    read <- character(length = read_length) #initialize read
    for (r in 1:read_length){ #simulate methylation of CpGs on read using reference parameters
      p <- rbeta(1, ref_sample["eta"], ref_sample["rho"])
      rand <- rbinom(1,1,p)
      read[r] <- ifelse(rand==1, "C", "T")
    }
    #output read_idx, cell_idx, marker_idx, read, num_obs
    output <- c(
      marker.index = as.numeric(ref_sample["marker_idx"]), 
      read = paste0(read, collapse = ""),
      nobs = num_obs,
      target = as.numeric(ref_sample["cell_idx"])
    )
    return(output)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(marker.index = as.numeric(marker.index)) %>%
    dplyr::mutate(nobs = as.numeric(nobs))
  
  
  ### EM Deconvolution
  weighted.psi.mat <- matrix(nrow = num_reads, ncol = num_celltypes)
  colnames(weighted.psi.mat) <- names(ref_list)
  psi.mat <-  matrix(nrow = num_reads, ncol = num_celltypes)
  colnames(psi.mat) <- names(ref_list)
  saved_weights <- list()
  for(i in 1:num_reads){
    # Calculate r
    r.vec = encode_binary(sim_pat$read[i])
    
    # Extract shape1, shape2, and beta.f for each marker that's associated with the ith read
    # Length of each of these vectors is the # of cell types
    shape1 = sapply(ref_list, function(x) return(x$shape1[x$marker.index==sim_pat$marker.index[i]]))
    shape2 = sapply(ref_list, function(x) return(x$shape2[x$marker.index==sim_pat$marker.index[i]]))
    beta.f = sapply(ref_list, function(x) return(x$beta.f[x$marker.index==sim_pat$marker.index[i]]))
    psi.vec = unlist(sapply(ref_list, function(x) return(x$psi.init[x$marker.index==sim_pat$marker.index[i]])))
    mu = sapply(ref_list, function(x) return(x$mu[x$marker.index==sim_pat$marker.index[i]]))
    sigma = sapply(ref_list, function(x) return(x$sigma[x$marker.index==sim_pat$marker.index[i]]))
    
    meth.frac = sum(r.vec) / length(r.vec) 
    weight.vec <- 1 / ceiling(abs(meth.frac - mu) / (num_error_sd*sigma))
    
    
    # For all cell types, iterate through each cell type c
    for(c in 1:length(psi.vec)){
      if(psi.vec[c]!=0){
        # For each CpG site r.i in read r (represented in r.vec)
        for(r.i in r.vec){
          # Compute P from the beta function
          P = beta(r.i + shape1[[c]], 1-r.i+shape2[[c]])/beta.f[[c]]
          # Update the psi value for each cell type by multiplying by P
          psi.vec[c] <- psi.vec[c] * P
        }
      }
    }
    # Output updated psi value into psi.mat
    weighted.psi.mat[i,] <- psi.vec * weight.vec
    psi.mat[i,] <- psi.vec #* weight.vec
    saved_weights[[i]] <- unique(weight.vec)
  }
  
  if(verbose) message("Performing EM")
  
  ### TODO: Loop and set different priors
  
  # initialize alpha prior
  #alpha <- rep(1/length(ref_list), length(ref_list))
  #alpha <- c(0.1-1e-7, 0.2, 0.3, 0.4, 1e-7)
  #alpha <- as.numeric(table(pat$V5) / sum(table(pat$V5)))
  num_of_inits <- 100
  alpha.inits <- lapply(1:num_of_inits, function(x){
    alpha <- runif(num_celltypes)
    alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
    return(alpha)
  })
  
  
  weighted.alphas <- lapply(alpha.inits, function(alpha){
    # Update alpha
    max.iter = 100
    i.iter = 1
    mad = vector()
    mad[1] <- 1
    all.alphas <- list()
    all.alphas[[i.iter]] <- alpha
    tol = 1e-3
    while(i.iter < max.iter & mad[i.iter] > tol){
      # Save old alpha
      alpha.old = alpha
      # Calculate new phi
      phi <- (alpha.old * weighted.psi.mat)/rowSums(alpha.old * weighted.psi.mat)
      # Calculate new alpha
      cs = colSums(phi * sim_pat$nobs)
      # cs = colSums(phi)
      new.alpha <- cs / sum(cs)
      
      i.iter = i.iter + 1
      
      # Check our threshold of mean absolute difference
      mad[i.iter] = mean(abs(alpha.old - new.alpha))/mean(new.alpha)
      alpha = new.alpha
      all.alphas[[i.iter]] <- alpha
    }
    return(alpha)
  })
  
  unweighted.alphas <- lapply(alpha.inits, function(alpha){
    # Update alpha
    max.iter = 100
    i.iter = 1
    mad = vector()
    mad[1] <- 1
    all.alphas <- list()
    all.alphas[[i.iter]] <- alpha
    tol = 1e-3
    while(i.iter < max.iter & mad[i.iter] > tol){
      # Save old alpha
      alpha.old = alpha
      # Calculate new phi
      phi <- (alpha.old * psi.mat)/rowSums(alpha.old * psi.mat)
      # Calculate new alpha
      cs = colSums(phi * sim_pat$nobs)
      # cs = colSums(phi)
      new.alpha <- cs / sum(cs)
      
      i.iter = i.iter + 1
      
      # Check our threshold of mean absolute difference
      mad[i.iter] = mean(abs(alpha.old - new.alpha))/mean(new.alpha)
      alpha = new.alpha
      all.alphas[[i.iter]] <- alpha
    }
    return(alpha)
  })
  
  return(c("weighted" = caret::RMSE(weighted.alphas[[1]], obs = true_alpha[,"alpha"]),
           "unweighted" = caret::RMSE(unweighted.alphas[[1]], obs = true_alpha[,"alpha"])))
}, cl = 6)
