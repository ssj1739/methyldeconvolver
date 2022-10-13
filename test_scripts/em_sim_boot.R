#library(methyldeconvolveR)
library(pbapply)
library(reshape2)
library(dbplyr)

#functions
encode_binary <- function(read){
  # Split string
  read.sp <- unlist(strsplit(read, split = ""))
  # Removing the missing values
  read.sp <- read.sp[read.sp!="."]
  
  # Encoding C as 1 and T as 0
  read.sp <- gsub(x = read.sp, "T", 0)
  read.sp <- gsub(x = read.sp, "C", 1)
  
  return(as.numeric(read.sp))
}

set.seed(1234)

###Simulation psuedo code
#input(num_celltypes, num_markers_per_celltype, num_reads)
num_celltypes = 5
num_markers_per_celltype = 50
num_reads = 10000


###Generate reference, simulated reads, and run full EM Deconvolution
###Output cell-type proportion estimates

#Generate simulated alpha:
alpha <- runif(num_celltypes)
alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
#hard code alphas:
#alpha <- c(0.75,0.15,0.03,0.02,0.05)
cell <- seq(1,num_celltypes)
true_alpha <- as.data.frame(cbind(cell, alpha))
colnames(true_alpha) <- c('cell_idx', 'alpha')
true_alpha$cell_idx <- as.factor(true_alpha$cell_idx)

noise = runif(1, 0.1, 0.67)
#Generate marker target indicator variable column:
sequence <- c((rep(0, num_celltypes-1)),1)
permutations <- arrangements::permutations(sort(unique(sequence)), freq = table(sequence),layout="list")
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
  ref[i, 'mu'] <- ta_ref['eta']/(ta_ref['eta']+ta_ref['rho'])
  ref[i, 'sigma'] <- sqrt((ta_ref['eta']*ta_ref['rho'])/(((ta_ref['eta']+ta_ref['rho'])^2)*(ta_ref['eta']+ta_ref['rho']+1)))
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
  } else {ref_sample <- ref_sub_noise[sample(nrow(ref_sub_noise), 1), ]}
  
  read_length <- round(runif(1,4,12)) #simulate read length
  num_obs <- rpois(1,10) #simulate number of observations of unique read
  read <- character(length = read_length) #initialize read
  for (r in 1:read_length){ #simulate methylation of CpGs on read using reference parameters
    p <- rbeta(1, ref_sample["eta"], ref_sample["rho"])
    rand <- rbinom(1,1,p)
    read[r] <- ifelse(rand==1, "C", "T")
  }
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

num_of_inits <- 100
alpha.inits <- lapply(1:num_of_inits, function(x){
  alpha <- runif(num_celltypes)
  alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
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

return(c("rmse" = caret::RMSE(unweighted.alphas[[1]], obs = true_alpha[,"alpha"]),
         "est_alpha"= list("unweighted.alphas" = unweighted.alphas[[1]])))

##Bootstrap SE & CI
n.boots <- 1000 #user input, should default to 10,000 but takes a long time

boot_output <- pblapply(1:n.boots, function(boot.sample){ # Loop through multiple bootstrap samples to compute SE
  #Sample boot.sample from sim_pat file  
  boot.sample <- sim_pat %>% slice_sample(n= num_reads, replace=TRUE)
    
    #Run Deconvolution on boot.sample
    ### EM Deconvolution
    psi.mat.boot <-  matrix(nrow = num_reads, ncol = num_celltypes)
    y.mat <- matrix(nrow = num_reads, ncol = 2, dimnames = dimnames(list(c(), c("y", "target"))))
    colnames(psi.mat.boot) <- names(ref_list)
    for(i in 1:num_reads){
      # Calculate r
      r.vec = encode_binary(boot.sample$read[i])
      
      # Extract shape1, shape2, and beta.f for each marker that's associated with the ith read
      # Length of each of these vectors is the # of cell types
      shape1 = sapply(ref_list, function(x) return(x$shape1[x$marker.index==boot.sample$marker.index[i]]))
      shape2 = sapply(ref_list, function(x) return(x$shape2[x$marker.index==boot.sample$marker.index[i]]))
      beta.f = sapply(ref_list, function(x) return(x$beta.f[x$marker.index==boot.sample$marker.index[i]]))
      psi.vec = unlist(sapply(ref_list, function(x) return(x$psi.init[x$marker.index==boot.sample$marker.index[i]])))
      mu = sapply(ref_list, function(x) return(x$mu[x$marker.index==boot.sample$marker.index[i]]))
      sigma = sapply(ref_list, function(x) return(x$sigma[x$marker.index==boot.sample$marker.index[i]]))
      target = sapply(ref_list, function(x) return(x$target[x$marker.index==boot.sample$marker.index[i]]))
      
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
      psi.mat.boot[i,] <- psi.vec 
    }
    
    #if(verbose) message("Performing EM")
    
    num_of_inits <- 10
    alpha.inits <- lapply(1:num_of_inits, function(x){
      alpha <- runif(num_celltypes)
      alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
      return(alpha)
    })
    
    boot.alphas <- pblapply(alpha.inits, function(alpha.i){
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
        transpose_prod.boot <- sweep(psi.mat.boot, MARGIN = 2, alpha.old, '*')
        phi.boot <- transpose_prod.boot/rowSums(transpose_prod.boot)
        # Calculate new alpha
        cs.boot = colSums(phi.boot * boot.sample$nobs)
        # cs = colSums(phi)
        new.alpha <- cs.boot / sum(cs.boot)
        
        i.iter = i.iter + 1
        
        # Check our threshold of mean absolute difference
        mad[i.iter] = mean(abs(alpha.old - new.alpha))/mean(new.alpha)
        alpha = new.alpha
        all.alphas[[i.iter]] <- alpha
      }
      return(alpha)
    })
    
    return("boot_alpha"= list("boot.alphas" = boot.alphas[[1]]))})
    

#Output cell type proportion estimates and create box plot
boot.out <- matrix(nrow = (n.boots), ncol = num_celltypes)
idx=1
for(i in 1:length(boot_output)){
    for(c in 1:num_celltypes){
      boot.out[idx,c] <- boot_output[[i]]$boot.alphas[c]}
    idx = idx+1
}
#compute bootstrap alpha estimate, se, 95% ci
alpha.ci <- 0.05 #can add user specified option for width of bootstrap ci, default to 95%
boot.out <- melt(boot.out)
colnames(boot.out) <- c('trial','cell_idx','alpha')
boot.summary <- boot.out %>% group_by(cell_idx)%>%summarise(mean = mean(alpha), sd = sd(alpha),se= sd/sqrt(n.boots),`2.5%`=quantile(alpha, 0.5*alpha.ci), `97.5%`=quantile(alpha,1-0.5*alpha.ci))
boot.summary$true_alpha <- true_alpha$alpha
