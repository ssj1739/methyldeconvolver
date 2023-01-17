# Simulate data

simulate_reference <- function(num_celltypes = 5,
                               num_markers_per_celltype = 25
                               ){
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
  
  return(ref_list)
  
  marker_df <- data.frame(
    chr = paste0("chr", sample(1:26, size = num_markers_per_celltype*num_celltypes, replace = T)),
    start = sample(1e4:2e6, size = num_markers_per_celltype*num_celltypes),
  )
}


simulate_pat <- function(ref_list, num_reads = 1e6, true_alpha = NA, zeta0 = -2, zeta1 = 5){
  if(length(true_alpha) == 1 | !is.data.frame(true_alpha)){
    message("True alpha not specified correctly. Simulating true alpha:")
    num_celltypes = length(ref_list)
    alpha <- runif(num_celltypes)
    alpha  <- unlist(lapply(alpha, function(x) x/sum(alpha) ))
    #hard code alphas:
    #alpha <- c(0.75,0.15,0.03,0.02,0.05)
    cell <- seq(1,num_celltypes)
    true_alpha <- as.data.frame(cbind(cell, alpha))
    colnames(true_alpha) <- c('cell_idx', 'alpha')
    true_alpha$cell_idx <- as.factor(true_alpha$cell_idx)
    print(true_alpha)
  }
  
  # flatten ref for compatibility
  ref = dplyr::bind_rows(ref_list, .id = "cell_idx") %>%
    mutate(cell_idx = as.numeric(cell_idx))
  
  #Generate reads for sample pat
  sim_pat <- pblapply(1:num_reads, function(x){
    cell_sim <- rmultinom(1,1,true_alpha[,"alpha"]) #simulate cell of origin
    cell_origin <- which(cell_sim == 1) 
    ref_sub_target <- subset(ref, cell_idx == cell_origin & target == 1) #subset reference by cell of origin
    ref_sub_noise <- subset(ref, cell_idx == cell_origin & target == 0) #subset reference by cell of origin
    noise <- exp(zeta0+zeta1*true_alpha[cell_origin,"alpha"])/(1+exp(zeta0+zeta1*true_alpha[cell_origin,"alpha"]))
    n <- rbinom(1,1,noise)
    if (n == 1){
      ref_sample <- ref_sub_target[sample(nrow(ref_sub_target), 1), ] #sample row from reference subset for alignment
    } else {ref_sample <- ref_sub_noise[sample(nrow(ref_sub_noise), 1), ]}
    
    read_length <- round(runif(1,4,12)) #simulate read length
    num_obs <- rpois(1,10) #simulate number of observations of unique read
    read <- character(length = read_length) #initialize read
    for (r in 1:read_length) { #simulate methylation of CpGs on read using reference parameters
      p <- rbeta(1, ref_sample[["shape1"]], ref_sample[["shape2"]])
      rand <- rbinom(1,1,p)
      read[r] <- ifelse(rand==1, "C", "T")
    }
    output <- data.frame(
      marker.index = as.numeric(ref_sample["marker.index"]), 
      read = paste0(read, collapse = ""),
      nobs = num_obs,
      target = as.numeric(ref_sample[["cell_idx"]])
    )
    return(output)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(marker.index = as.numeric(marker.index)) %>%
    dplyr::mutate(nobs = as.numeric(nobs))
  

  
  return(sim_pat)
}
