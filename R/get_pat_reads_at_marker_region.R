#' get_pat_reads_at_marker_region
#'
#' @param marker 
#' @param pat 
#' @param pat.dir 
#' @param cell_type 
#'
#' @return
#' @export
#'
get_pat_reads_at_marker_region <- function(marker, marker.ind, pat.file=NULL, pat.dir=NULL, cell_type=NULL, 
                                           filter.noninf = T, filter.length = 3, filter.inf.length = 3){
  # TODO: Get a vector of the pat reads from a cell_type at a specific marker region.
  # First, read files:
  if(!is.null(pat.dir)){
    if(is.null(cell_type)){
      stop("Can't specify pat.dir without also specifying cell_type.")
    }
    
    pat.files = dir(pattern = "*.pat.gz$", pat.dir, full.names = F)
    pat.cell_types <- tolower(sapply(pat.files, function(x) strsplit(x, split = "_")[[1]][1]))
    if(!cell_type %in% pat.cell_types){
      stop("Specified cell_type is not found in the pat.dir.")
    }
    # Read PAT files for specified cell type
    pc_pat.files <- pat.files[pat.cell_types == cell_type]
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
    
    final_pat <- dplyr::bind_rows(pc_pat.list)
  }else{
    if(!is.null(pat.file)){
      pat <- read_pat(path = pat.file, verbose = T)
      final_pat <- pat
    }else{
      stop("Neither pat.file nor pat.dir was specified.")
    }
  }
  
  # Second, compute overlap
  overlap <- overlap_marker_pat(pat = final_pat, marker = marker, split_reads = F)
  
  # Third, pull pat reads from specific marker.ind
  which.pat <- overlap$overlaps@from[overlap$overlaps@to == marker.ind]
  if(!length(which.pat)>0){
    warning("No PAT reads that align to this marker region.")
    return(NULL)
  }
  
  reads <- rep(final_pat$read[which.pat], times = final_pat$nobs[which.pat])
  return(reads)
}
