# Workflow for reference

#' learn_reference
#' @description Learn a reference set from given reference pat files of cell types of interest.
#'
#' @param marker.file Path to marker file
#' @param pat.dir Path to directory containing reference pat files
#' @param save.output character - desired path to output reference, saved as serialized R object (.rds)
#'
#' @return reference object
#' @export
#'
#' @examples
#' \dontrun{
#' learn_reference(marker.file = "marker.txt", pat.dir = "data/ref/", save.output = "reference.rds")
#' }
learn_reference <- function(marker.file, pat.dir, save.output = "", verbose = F){
  
  # Preprocess marker, ensure correct format
  if(verbose) message("Reading in marker file")
  marker = read_marker(marker.file)
  
  # Identify # of Pat files in pat.dir, ensure correct formats
  pat.files = dir(pattern = "*.pat*", pat.dir, full.names = F)
  
  # Loop through each pat file in directory
  if(verbose) message("Starting to loop through pat files:")
  beta_celltype_fits <- list()
  pat.num <- 1
  for(pf in pat.files){
    if(verbose) message(paste0("Reading ", pat.num, " of ", length(pat.files), " PAT files"))
    pat <- read_pat(path = paste0(pat.dir,"/",pf), verbose = verbose)
    pat.num <- pat.num+1
    
    if(verbose) message("Checking overlap of ", pf)
    overlap <- overlap_marker_pat(pat = pat, marker = marker)
    
    if(verbose) message("Fitting beta distributions.")
    pf_name <- tolower(strsplit(pf, split = ".pat", fixed = T)[[1]][1])
    beta_celltype_fits[[pf_name]] <- fit_beta(overlaps.list = overlap)
  }
  
  # Return reference format, with markers and shape params/beta.f for each celltype
  marker.subset <- marker %>%
    select(chr, startCpG, endCpG, target)
  
  bad.markers <- c()
  for(cell_type in names(beta_celltype_fits)){
  
    beta_celltype_fits_subset <- beta_celltype_fits[[cell_type]][marker.subset$target==cell_type,]
    bad.markers <- c(bad.markers, beta_celltype_fits_subset$marker.index[beta_celltype_fits_subset$psi.init==0])
  }
  bad.markers <- unique(bad.markers)
  
  
  output <- list(marker = marker, beta_celltype_fits = beta_celltype_fits)
  
  if(save.output!=""){
    if(is.logical(save.output) & isTRUE(save.output))
      save.output <- paste0("deconv_reference_",Sys.Date(),".rds")
    saveRDS(object = output, file = save.output)
  }
  
  return(output)
}
