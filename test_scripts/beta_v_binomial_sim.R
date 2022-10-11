library(pbapply)
library(reshape2)

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
num_celltypes = 20
num_markers_per_celltype = 50
num_reads = 10000


#Generate simulated alpha:
alpha <- runif(num_celltypes)
alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
#hard code alphas:
#alpha <- c(0.75,0.15,0.03,0.02,0.05)
cell <- seq(1,num_celltypes)
true_alpha <- as.data.frame(cbind(cell, alpha))
colnames(true_alpha) <- c('cell_idx', 'alpha')
true_alpha$cell_idx <- as.factor(true_alpha$cell_idx)

n.iters <- 100
n.trials <- 3

trial_output <- pblapply(1:n.trials, function(trial){ # Loop through multiple trials with varying amounts of noise
  trial.iters <- pblapply(1:n.iters, function(x){ #perform multiple iterations for each noise level
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
    {
      if(trial == 1){
      for (i in 1:nrow(ref)){
        if (ref[i,'target'] == 1){
          ref[i,'eta'] <- runif(1, min = 0.001, max = 2)
          ref[i, 'rho'] <- runif(1,8,20)
          } else {
          ref[i,'eta'] <- runif(1,8,20)
          ref[i, 'rho'] <- runif(1,0.001,2)
        }
        }
      }else if(trial == 2){
        for (i in 1:nrow(ref)){
          if (ref[i,'target'] == 1){
            ref[i,'eta'] <- runif(1, min = 1.5, max = 3)
            ref[i, 'rho'] <- runif(1,6,10)
            }else{
              ref[i,'eta'] <- runif(1,6,10)
              ref[i, 'rho'] <- runif(1,1.5,3)
              }
          }
        } else {
          for (i in 1:nrow(ref)){
            if (ref[i,'target'] == 1){
              ref[i,'eta'] <- runif(1, min = 2, max = 8)
              ref[i, 'rho'] <- runif(1,2,8)
              } else {
                ref[i,'eta'] <- runif(1,8,20)
                ref[i, 'rho'] <- runif(1,0.001,2)
              }
          }
        }
    }
    
    #Calculate mean and SD of methylation pattern at marker target celltype (in ref)
    # aka average methylation for all markers of a given cell type, and save in ref
      for(i in 1:nrow(ref)){
        ref[i, 'mu'] <- ref[i,'eta']/(ref[i,'eta']+ref[i,'rho'])
        ref[i, 'sigma'] <- sqrt((ref[i,'eta']*ref[i,'rho'])/(((ref[i,'eta']+ref[i,'rho'])^2)*(ref[i,'eta']+ref[i,'rho']+1)))
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
    #y.mat <- matrix(nrow = num_reads, ncol = 2, dimnames = dimnames(list(c(), c("y", "target"))))
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
      #meth.frac = sum(r.vec) / length(r.vec) 
      # if meth.frac < mu+(2*sigma) & meth.frac > mu-(2*sigma) then y.vec=1, else y.vec=0 
      #n.sigma = 3
      #upper_bound <- mu[which(target==1)]+n.sigma*sigma[which(target==1)]
      #lower_bound <-  mu[which(target==1)]-n.sigma*sigma[which(target==1)]
      #y.mat[i,1] <- as.numeric(meth.frac < upper_bound & meth.frac > lower_bound)
      #y.mat[i,2] <- which(target==1)
      
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
    
    psi.mat.beta <-  matrix(nrow = num_reads, ncol = num_celltypes)
    #y.mat <- matrix(nrow = num_reads, ncol = 2, dimnames = dimnames(list(c(), c("y", "target"))))
    colnames(psi.mat.beta) <- names(ref_list)
    for(i in 1:num_reads){
      # Calculate r
      r.vec = encode_binary(sim_pat$read[i])
      
      # Extract shape1, shape2, and beta.f for each marker that's associated with the ith read
      # Length of each of these vectors is the # of cell types
      shape1 = sapply(ref_list, function(x) return(x$shape1[x$marker.index==sim_pat$marker.index[i]]))
      shape2 = sapply(ref_list, function(x) return(x$shape2[x$marker.index==sim_pat$marker.index[i]]))
      beta.f = sapply(ref_list, function(x) return(x$beta.f[x$marker.index==sim_pat$marker.index[i]]))
      psi.vec.beta = unlist(sapply(ref_list, function(x) return(x$psi.init[x$marker.index==sim_pat$marker.index[i]])))
      mu = sapply(ref_list, function(x) return(x$mu[x$marker.index==sim_pat$marker.index[i]]))
      sigma = sapply(ref_list, function(x) return(x$sigma[x$marker.index==sim_pat$marker.index[i]]))
      target = sapply(ref_list, function(x) return(x$target[x$marker.index==sim_pat$marker.index[i]]))
      
      #Binary is read consistent with marker target?
      #meth.frac = sum(r.vec) / length(r.vec) 
      # if meth.frac < mu+(2*sigma) & meth.frac > mu-(2*sigma) then y.vec=1, else y.vec=0 
      #n.sigma = 3
      #upper_bound <- mu[which(target==1)]+n.sigma*sigma[which(target==1)]
      #lower_bound <-  mu[which(target==1)]-n.sigma*sigma[which(target==1)]
      #y.mat[i,1] <- as.numeric(meth.frac < upper_bound & meth.frac > lower_bound)
      #y.mat[i,2] <- which(target==1)
      
      # For all cell types, iterate through each cell type c
      for(c in 1:length(psi.vec.beta)){
        if(psi.vec.beta[c]!=0){
          # sample CpG site r.i in read r (represented in r.vec)
            r.i <- sample(r.vec, 1, replace=FALSE)
            # Compute P from the binomial function
            P = dbinom(r.i,1,mu[[c]])
            # Set the psi value for each cell type to P
            psi.vec.beta[c] <-  P
        }
      }
      # Output updated psi value into psi.mat
      psi.mat.beta[i,] <- psi.vec.beta 
    }
    
    #if(verbose) message("Performing EM")
    
    num_of_inits <- 1
    alpha.inits <- lapply(1:num_of_inits, function(x){
      alpha <- runif(num_celltypes)
      alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
      return(alpha)
    })
    
    
    beta.alphas <- pblapply(alpha.inits, function(alpha.i){
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
        transpose_prod <- sweep(psi.mat.beta, MARGIN = 2, alpha.old, '*')
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
    
    return(c("trial"=trial,"beta.rmse" = caret::RMSE(beta.alphas[[1]], obs = true_alpha[,"alpha"]),
             "unweighted.rmse" = caret::RMSE(unweighted.alphas[[1]], obs = true_alpha[,"alpha"]),
             list("beta.alphas" = beta.alphas[[1]]), 
             list("unweighted.alphas" = unweighted.alphas[[1]])))
  })
  return(list(trial.iters))
}, cl = 6)

#Output cell type proportion estimates and create box plot
beta.sim.out <- matrix(nrow = (n.iters*n.trials), ncol = num_celltypes)
ref.trial.vec <- c()
rmse.vec <- c()
idx=1
for(i in 1:length(trial_output)){
  for(n in 1:n.iters){
    for(c in 1:num_celltypes){
      beta.sim.out[idx,c] <- trial_output[[i]][[1]][[n]]$unweighted.alphas[c]}
    ref.trial.vec <-append(ref.trial.vec,trial_output[[i]][[1]][[n]]$trial)
    rmse.vec <-append(rmse.vec,trial_output[[i]][[1]][[n]]$unweighted.rmse)
    idx = idx+1
  }
}
beta.sim.out <- melt(beta.sim.out)
colnames(beta.sim.out) <- c('trial','cell_idx','alpha')
beta.sim.out$model <- rep('beta',length(beta.sim.out$alpha))
beta.sim.out$rmse <- rep(rmse.vec, num_celltypes)

binomial.sim.out <- matrix(nrow = (n.iters*n.trials), ncol = num_celltypes)
idx=1
rmse.vec <- c()
for(i in 1:length(trial_output)){
  for(n in 1:n.iters){
    for(c in 1:num_celltypes){
      binomial.sim.out[idx,c] <- trial_output[[i]][[1]][[n]]$beta.alphas[c]}
    rmse.vec <-append(rmse.vec,trial_output[[i]][[1]][[n]]$beta.rmse)
    idx = idx+1
  }
}

binomial.sim.out <- melt(binomial.sim.out)
colnames(binomial.sim.out) <- c('trial','cell_idx','alpha')
binomial.sim.out$model <- rep('binomial',length(binomial.sim.out$alpha))
binomial.sim.out$rmse <- rep(rmse.vec, num_celltypes)

sim.out <- rbind(beta.sim.out,binomial.sim.out)
sim.out$model <- as.factor(sim.out$model)
sim.out$cell_idx <- as.factor(sim.out$cell_idx)
sim.out$ref.trial <- rep(ref.trial.vec, (num_celltypes*2))

sim.out <- within(sim.out, {   
  marker_shape <- NA # need to initialize variable
  marker_shape[ref.trial == 1] <- "Skewed"
  marker_shape[ref.trial == 2] <- "Hemi"
  marker_shape[ref.trial == 3] <- "Mixed"
} )

write.csv(sim.out,"C:/Users/Patri/Documents/School/Georgetown/Thesis/Deconv Sim/sim.out.20.50.beta_v_binomial.csv", row.names = FALSE)

p <- ggplot(data = sim.out) + 
  geom_boxplot(aes(x=cell_idx, y=alpha, group=interaction(cell_idx,marker_shape), fill=marker_shape, col=marker_shape),alpha=0.3) +
  labs(title="20 Cell Types - 50 Markers Per - Beta vs Binomial",
       x ="Cell Type", y = "Proportion Estimate", col="Ref Marker Shape")+
  guides(size = "legend", fill = "none")+
  facet_wrap(~model)+
  theme_minimal() +
  #theme(legend.position="none")+
  geom_point(data=data.frame(true_alpha), aes(x=cell_idx, y=alpha), color="black", size=2)
p

summary <- sim.out %>% group_by(model,noise.cat,cell_idx)%>%summarise(mean = mean(alpha), sd = sd(alpha))
rmse <- sim.out %>% group_by(model,noise.cat)%>%summarise(rmse = mean(rmse))
