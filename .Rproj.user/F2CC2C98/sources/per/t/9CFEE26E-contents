#'@description Merging the PAT and Marker files together
#'@param pat Some description of PAT files
#'@param marker Some description of Marker files
#'
#'@export
#'
overlap_marker_pat <- function(pat, marker, minoverlap = 1, type = "any", ...){
  require(GenomicRanges)
  require(dplyr)
  
  # Convert marker to GRanges object:
  marker.ranges <- GenomicRanges::makeGRangesFromDataFrame(marker %>% dplyr::select(-c("start", "end")), 
                                                           start.field = "startCpG", 
                                                           end.field = "endCpG", 
                                                           keep.extra.columns = T)
  
  # Convert pat to GRanges object:
  pat$end <- pat$start + nchar(pat$read) # Add column for end to load into GR
  pat.ranges <- GenomicRanges::makeGRangesFromDataFrame(pat, keep.extra.columns = T)
  
  # Find overlaps, taking into account overlap
  marker.pat.overlaps <- findOverlaps(query = pat.ranges, subject = marker.ranges, minoverlap = 1, type = "any", )

  return(list(overlaps = marker.pat.overlaps, pat.gr = pat.ranges, marker.gr = marker.ranges))
}
