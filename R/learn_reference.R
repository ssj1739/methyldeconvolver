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
learn_reference <- function(marker.file, 
                            pat.dir, 
                            save.output = "", 
                            verbose = T, 
                            n_threads = 1,
                            split_reads = F,
                            filter.inf.length = 3,
                            filter.length = 3,
                            filter.noninf = T,
                            min.cell.param = NA, ...){
  
  # Preprocess marker, ensure correct format
  if(verbose) message("Reading in marker file")
  marker = read_marker(marker.file, ...)
  
  # Identify # of Pat files in pat.dir, ensure correct formats
  pat.files = dir(pattern = "*.pat.gz$", pat.dir, full.names = F)
  pat.cell_types <- sapply(pat.files, function(x) strsplit(x, split = "_")[[1]][1]) %>%
    tolower() %>% 
    gsub(pattern = "[[:punct:]]", replacement = "") %>%
    gsub(pattern = " ", replacement = "")
  
  # Verify that all PAT file cell types are contained in the marker file "targets" column
  missing_pats <- which(!tolower(pat.cell_types) %in% tolower(marker$label)) # Use `tolower` to normalize case mismatch
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
  if(length(pat.files) == 0)
    stop("ERROR: No PAT files found in specified pat.dir.")
  beta_celltype_fits <- list()
  
  cell_type.num <- 1
  
  for(pc in unique(pat.cell_types)){
    if(verbose) message(paste0("Reading PAT files from cell-type ",cell_type.num," out of ", length(unique(pat.cell_types)), ": ", pc))
    pc_pat.files <- pat.files[pat.cell_types == pc]
    pat.num <- 1
    
    pc_pat.list <- list()
    #overlap_list <- list()
    for(pf in pc_pat.files){
      if(verbose) message(paste0("Reading ", pat.num, " of ", length(pc_pat.files), " PAT files"))
      pat <- read_pat(path = paste0(pat.dir,"/",pf),
                      verbose = verbose,
                      filter.noninf = filter.noninf,
                      filter.length = filter.length,
                      filter.inf.length = filter.inf.length) # Filter out reads in reference PAT containing less than 3 CpGs

      pat.num <- pat.num+1
      pc_pat.list[[pf]] <- pat
    }
    
    pc_pat.merged <- dplyr::bind_rows(pc_pat.list)
    #rm(pc_pat.list)

    if(verbose) message("Checking overlap of PAT files from ", pc)
    overlap <- overlap_marker_pat(pat = pc_pat.merged,
                                  marker = marker,
                                  n_threads = n_threads,
                                  split_reads = split_reads)

    #Learn beta distribution for each maker region
    if(verbose) message("Fitting beta distributions.")
    beta_celltype_fits[[pc]] <- fit_beta(overlaps.list = overlap)
    cell_type.num <- cell_type.num + 1
    rm(pc_pat.list, pc_pat.merged, overlap)
    gc(full = T)
  }
  
  # Return reference format, with markers and shape params/beta.f for each celltype
   marker.subset <- marker
  #   as.data.frame() %>%
  #   select(chr, startCpG, endCpG, target)
  
  #Identify cell/marker pairs with insufficient reference data or unable to fit beta     
  #TODO: Make sure this works!
  if(is.na(min.cell.param)){
    # by default, set to 20% of cell types
    min.cell.param <- round(0.2 * length(beta_celltype_fits))
  }
   
  psi.mat.init <- sapply(beta_celltype_fits, function(x) return(x$psi.init)) 
  bad.markers <- which(rowSums(psi.mat.init) < min.cell.param) # markers where no information was learned for any cell type
  marker.out <- marker[-bad.markers,]
  beta_celltype_fits.out <- lapply(beta_celltype_fits, function(x) x[-bad.markers,])
  
  output <- list(marker = marker.out, beta_celltype_fits = beta_celltype_fits.out)
  
  if(save.output!=""){
    if(is.logical(save.output) & isTRUE(save.output))
      save.output <- paste0("deconv_reference_",Sys.Date(),".rds")
    saveRDS(object = output, file = save.output)
  }
  
  return(output)
}
