# Workflow for reference

#' learn_reference
#' @description Learn a reference set from given reference pat files of cell types of interest.
#'
#' @param marker.file Path to marker file
#' @param pat.dir Path to directory containing reference pat files. See notes.
#' @param save.output character - desired path to output reference, saved as serialized R object (.rds)
#'
#' @return reference object
#' @export
#'
#' @note Name formats of the reference PAT files should follow the following convention:
#' - Files should end in .pat.gz (bgzipped PAT files outputted from wgbstools).
#' - Files should have the name of the cell type (which perfectly matches the name in the marker file) 
#' followed by an underscore, followed by any other names.
#' e.g. Bcell_1000_ExampleCellGeneratedInHouse.pat.gz
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
  pat.cell_types <- tolower(sapply(pat.files, function(x) strsplit(x, split = "_")[[1]][1]))
  
  # Verify that all PAT file cell types are contained in the marker file "targets" column
  missing_pats <- which(!tolower(pat.cell_types) %in% tolower(marker$target)) # Use `tolower` to normalize case mismatch
  if(length(missing_pats)>0){
    if(verbose){
      message("The following PAT files represent cell types not found in the marker file:")
      cat(pat.files[missing_pats], sep = '\n')
    }
    # Remove PAT files without a matching target cell type
    pat.files <- pat.files[-missing_pats]
    pat.cell_types <- pat.cell_types[-missing_pats]
  }
  
  # Loop through each pat file in directory
  if(verbose) message("Starting to loop through pat files:")
  beta_celltype_fits <- list()
  
  cell_type.num <- 1
  for(pc in unique(pat.cell_types)){
    if(verbose) message(paste0("Reading PAT files from cell-type: ", pc, "; ", cell_type.num, " out of ", length(unique(pat.cell_types))))
    pc_pat.files <- pat.files[pat.cell_types == pc]
    pc_pat.list <- list()
    pat.num <- 1
    for(pf in pc_pat.files){
      if(verbose) message(paste0("Reading ", pat.num, " of ", length(pc_pat.files), " PAT files"))
      pat <- read_pat(path = paste0(pat.dir,"/",pf), 
                      verbose = verbose, 
                      filter.noninf = T, 
                      filter.length = 3) # Filter out reads in reference PAT containing less than 3 CpGs
      pat.num <- pat.num+1
      pc_pat.list[[pf]] <- pat
    }
    
    pc_pat.merged <- bind_rows(pc_pat.list)
      
    if(verbose) message("Checking overlap of PAT files from ", pc)
    overlap <- overlap_marker_pat(pat = pc_pat.merged, marker = marker)
      
    if(verbose) message("Fitting beta distributions.")
    beta_celltype_fits[[pc]] <- fit_beta(overlaps.list = overlap)
    cell_type.num <- cell_type.num + 1
    rm(pc_pat.merged, overlap)
    gc(full = T)
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
