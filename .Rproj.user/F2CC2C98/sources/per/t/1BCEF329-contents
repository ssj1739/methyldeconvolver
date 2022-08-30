# Simulate reference

#' simulate_reference
#' @description Simulate a deconvolution reference.
#' Requires a number of cell types and number of markers per cell type
#' Creates a reference by simulating a methylation level for each marker
#' and indicating a target for each marker.
#'
#' @param num_celltypes 
#' @param num_markers_per_celltype 
#'
#' @return ref_list object (similar to learn_reference)
#' @export
#'
#' @examples
simulate_reference <- function(num_celltypes, num_markers_per_celltype){
  
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
  
  return(ref_list)
}
