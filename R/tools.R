#' encode_binary Calculate binary representation of methylation after bisulfite conversion
#'
#' @param read string of DNA methylation read (containing T and C)
#'
#' @description Given a string of T and C, calculates the binary representation of that methylation.
#' 0 represents non-methylated cytosines, and 1 represents methylated cytosines (converted to T)
#' @return integer vector
#' @export
#'
#' @examples
#' \dontrun{
#' encode_binary("C...CT")
#' }
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

#' filter_pat
#'
#' @param pat 
#' @param filter.noninf 
#' @param filter.length 
#' @param filter.inf.length 
#' @param output.as.granges 
#' @param verbose 
#'
#' @return
#' @export
#'
filter_pat <- function(pat, filter.noninf = T, 
                       filter.length = 3, 
                       filter.inf.length = 3,
                       output.as.granges = F, verbose = F){
  
  if("GRanges" %in% class(pat)){
    output.as.granges = T
    pat <- as.data.frame(pat) %>%
      dplyr::select(-strand)
  }
  
  # Filtering pat files by number of CpGs
  pat.filt <- pat %>%
    dplyr::filter(nchar(read) >= filter.length)
  
  # Filtering PAT files by read information (should contain CpG with known methylation status)
  if(isTRUE(filter.noninf)){
    pat.filt <- pat.filt %>%
      dplyr::filter(grepl("C|T", read))
  }
  
  # Filtering PAT files by number of informative CpGs
  pat.filt <- pat.filt %>%
    dplyr::filter(str_count(string = read, pattern = "C|T") > filter.inf.length)
  
  if(verbose){
    message("Finished filtering.")
    message(paste0("Filtered out ", sum(pat$nobs) - sum(pat.filt$nobs),
                   " (",
                   round(x = 100*(sum(pat$nobs) - sum(pat.filt$nobs))/sum(pat$nobs),digits = 2),
                   "%) observed reads out of ", sum(pat$nobs), "."))
  }
  
  if(output.as.granges){
    pat.filt <- pat.filt %>%
      dplyr::mutate(end = start + stringr::str_length(pat.filt$read)) # Add column for end to load into GR
    pat.ranges <- GenomicRanges::makeGRangesFromDataFrame(pat.filt, 
                                                          keep.extra.columns = T, 
                                                          ignore.strand = T, 
                                                          end.field = 'end', 
                                                          start.field = 'start')
    return(pat.ranges)
  }
  
  return(pat.filt)
}
