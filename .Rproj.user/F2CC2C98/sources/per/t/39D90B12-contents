library(methyldeconvolveR)
library(pbapply)
library(tidyverse)

###Simulation psuedo code
#input(num_celltypes, num_markers_per_celltype, num_reads)
# num_celltypes = 10
num_markers_per_celltype = 25
num_reads = 50000
num_celltypes = 5
#noise = 0.5

#Generate simulated alpha:
# alpha <- runif(num_celltypes)
# alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
# cell <- seq(1,num_celltypes)
# true_alpha <- cbind(cell, alpha)
# colnames(true_alpha) <- c('cell_idx', 'alpha')

trial_output <- pblapply(1:50, function(trial){ # Loop through multiple trials with varying amounts of noise
  # Trial for different levels of noise
  noise = 0.02*trial
  
  # Trial for different number of cell types

  alpha <- runif(num_celltypes)
  alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
  cell <- seq(1,num_celltypes)
  true_alpha <- cbind(cell, alpha)
  colnames(true_alpha) <- c('cell_idx', 'alpha')
  
  #Generate marker target indicator variable column:
  sequence <- c((rep(0, num_celltypes-1)),1)
  #permutations <- unique(arrangements::permutations(sequence, layout = "list"))
  permutations <- lapply(1:num_celltypes, function(x){
    perm <- rep(0, num_celltypes)
    perm[x] <- 1
    return(perm)
  })
  target <- rep(unlist(permutations), num_markers_per_celltype)
  
  #Simulate reference parameters, initialize columns
  total_markers <- num_celltypes*num_markers_per_celltype
  marker_idx <- rep(1:total_markers, each= num_celltypes)
  cell_idx <- rep(1:num_celltypes, total_markers)
  eta <- vector(mode = "numeric",total_markers*num_celltypes)
  rho <- vector(mode = "numeric", total_markers*num_celltypes)
  mu <- vector(mode = "numeric", total_markers*num_celltypes)
  sigma <- vector(mode = "numeric", total_markers*num_celltypes)
  ref <- cbind(marker_idx, cell_idx, target, eta, rho, mu, sigma)
  colnames(ref) <- c('marker_idx', 'cell_idx', 'target', 'eta', 'rho', 'mu', 'sigma')
  
  # Populate eta and rho with runif values (random beta distributions for each marker)
  for (i in 1:nrow(ref)){
    if (ref[i,'target'] == 1){
      ref[i,'eta'] <- runif(1, min = 0.001, max = 2)
      ref[i, 'rho'] <- runif(1,8,20)
    } else {
      ref[i,'eta'] <- runif(1,8,20)
      ref[i, 'rho'] <- runif(1,0.001,2)
    }
  }
  
  #Calculate mean and SD of methylation pattern at marker target celltype (in ref)
  # aka average methylation for all markers of a given cell type, and save in ref
  for(i in 1:nrow(ref)){
    mi <- ref[i,"marker_idx"]
    ta_ref <- ref[ref[,"marker_idx"]==mi & ref[,"target"]==1,]
    if(is.vector(ta_ref)){
      ref[i, 'mu'] <- ta_ref['eta']/(ta_ref['eta']+ta_ref['rho'])
      ref[i, 'sigma'] <- sqrt((ta_ref['eta']*ta_ref['rho'])/(((ta_ref['eta']+ta_ref['rho'])^2)*(ta_ref['eta']+ta_ref['rho']+1)))
    }else{
      ref[i, 'mu'] <- ta_ref[,'eta']/(ta_ref[,'eta']+ta_ref[,'rho'])
      ref[i, 'sigma'] <- sqrt((ta_ref[,'eta']*ta_ref[,'rho'])/(((ta_ref[,'eta']+ta_ref[,'rho'])^2)*(ta_ref[,'eta']+ta_ref[,'rho']+1)))
    }
    
  }
  
  # Transform ref into format used in package
  ref_list <- ref %>%
    as.data.frame() %>%
    dplyr::select(target=target, shape1=eta, shape2=rho, marker.index = marker_idx, mu = mu, sigma = sigma) %>%
    dplyr::mutate(beta.f = beta(shape1, shape2), psi.init = 1) %>%
    split(cell_idx)
  
  #Generate reads for sample pat
  sim_pat <- lapply(1:num_reads, function(x){
    cell_sim <- rmultinom(1,1,true_alpha[,"alpha"]) #simulate cell of origin
    cell_origin <- which(cell_sim == 1) 
    ref_sub_target <- subset(ref, cell_idx == cell_origin & target == 1) #subset reference by cell of origin
    ref_sub_noise <- subset(ref, cell_idx == cell_origin & target == 0) #subset reference by cell of origin
    n <- rbinom(1,1,1-noise)
    if (n == 1){
      ref_sample <- ref_sub_target[sample(nrow(ref_sub_target), 1), ] #sample row from reference subset for alignment
    } else {
      ref_sample <- ref_sub_noise[sample(nrow(ref_sub_noise), 1), ]
      }
    
    read_length <- round(runif(1,4,12)) #simulate read length
    num_obs <- rpois(1,10) #simulate number of observations of unique read
    read <- character(length = read_length) #initialize read
    for (r in 1:read_length){ #simulate methylation of CpGs on read using reference parameters
      p <- rbeta(1, ref_sample["eta"], ref_sample["rho"])
      rand <- rbinom(1,1,p)
      read[r] <- ifelse(rand==1, "C", "T")
    }
    output <- data.frame(
      marker.index = as.numeric(ref_sample["marker_idx"]), 
      read = paste0(read, collapse = ""),
      nobs = num_obs,
      target = as.numeric(ref_sample["cell_idx"])
    )
    return(output)
  }) %>%
    dplyr::bind_rows()
    
  ### EM Deconvolution
  psi.mat <-  matrix(nrow = num_reads, ncol = num_celltypes)
  y.mat <- matrix(nrow = num_reads, ncol = 2, dimnames = dimnames(list(c(), c("y", "target"))))
  colnames(psi.mat) <- names(ref_list)
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
    target = sapply(ref_list, function(x) return(x$target[x$marker.index==sim_pat$marker.index[i]]))
    
    #Binary is read consistent with marker target?
    meth.frac = sum(r.vec) / length(r.vec) 
    # if meth.frac < mu+(2*sigma) & meth.frac > mu-(2*sigma) then y.vec=1, else y.vec=0 
    n.sigma = 3
    upper_bound <- mu[which(target==1)]+n.sigma*sigma[which(target==1)]
    lower_bound <-  mu[which(target==1)]-n.sigma*sigma[which(target==1)]
    y.mat[i,1] <- as.numeric(meth.frac < upper_bound & meth.frac > lower_bound)
    y.mat[i,2] <- which(target==1)
    
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
    psi.mat[i,] <- psi.vec 
  }
  
  #if(verbose) message("Performing EM")
  
  num_of_inits <- 10
  alpha.inits <- lapply(1:num_of_inits, function(x){
    alpha <- runif(num_celltypes)
    alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
    return(alpha)
  })
  
  
  weighted.alphas <- pblapply(alpha.inits, function(alpha.i){
    # Update alpha
    max.iter = 100
    i.iter = 1
    mad = vector()
    mad[1] <- 1
    all.alphas <- list()
    all.alphas[[i.iter]] <- alpha.i
    tol = 1e-3
    
    while(i.iter < max.iter & mad[i.iter] > tol){
      # Save old alpha
      alpha.old = alpha
      #calculate zeta
      #fit logistic regression to y.vec given current estimate of alpha
      learn.zeta.mat <- data.frame(y = y.mat[,1], alpha = alpha[y.mat[,2]])
      zeta.fit <- glm(y ~ alpha, data = learn.zeta.mat, family="binomial")
        
      omega.vec <- predict(zeta.fit, newdata = data.frame(alpha = alpha), type = "response")
      
      # Calculate new phi
      transpose_prod <- sweep(psi.mat, MARGIN = 2, alpha.old * omega.vec, '*')
      #transpose_prod <- t((alpha.old * omega.vec) * t(psi.mat))
      phi <- transpose_prod/rowSums(transpose_prod)
      #phi.uw <- (alpha.old * psi.mat) / rowSums(alpha.old * psi.mat)
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
  
  unweighted.alphas <- pblapply(alpha.inits, function(alpha.i){
    # Update alpha
    max.iter = 100
    i.iter = 1
    mad = vector()
    mad[1] <- 1
    all.alphas <- list()
    all.alphas[[i.iter]] <- alpha.i
    tol = 1e-3
    while(i.iter < max.iter & mad[i.iter] > tol){
      # Save old alpha
      alpha.old = alpha
      # Calculate new phi
      transpose_prod <- sweep(psi.mat, MARGIN = 2, alpha.old, '*')

      phi <- transpose_prod/rowSums(transpose_prod)
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
  
  output_df <- data.frame(
    trial_num = trial,
    weighted = caret::RMSE(weighted.alphas[[1]], obs = true_alpha[,"alpha"]),
    unweighted = caret::RMSE(unweighted.alphas[[1]], obs = true_alpha[,"alpha"]))
  
  return(output_df)
}, cl = 6)

# Visualize trial output

trial_output_df <- trial_output %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(noise_rate = trial_num * 0.01) %>%
 # dplyr::mutate(num_celltypes = trial_num*2) %>%
  pivot_longer(cols = c(weighted, unweighted), names_to = "model_type", values_to = "RMSE")

ggplot(data = trial_output_df, aes(x = noise_rate, y = RMSE, color = model_type)) +
  geom_point()

ggsave("RMSE_noise.png", width = 7, height = 5)
