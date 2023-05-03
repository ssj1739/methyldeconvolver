
#' Overlap marker and PAT files
#' 
#'@description Merging the PAT and Marker files together
#'@param pat Some description of PAT files
#'@param marker Some description of Marker files
#'@param n_threads Number of threads to use (default = 1)
#'
#'@export
#'
overlap_marker_pat <- function(pat, marker, split_reads = F, n_threads = 1){
  require(GenomicRanges)
  require(dplyr)
  require(data.table)
  require(pbapply)
  
  # Convert marker to GRanges object:
  marker.ranges <- marker
  # marker.ranges <- GenomicRanges::makeGRangesFromDataFrame(marker %>% 
  #                                                            dplyr::select(-c("start", "end")), 
  #                                                          start.field = "startCpG", 
  #                                                          end.field = "endCpG", 
  #                                                          keep.extra.columns = T)
  
  # Convert pat to GRanges object:
  pat <- pat %>%
    mutate(end = start + stringr::str_length(pat$read)) # Add column for end to load into GR
  pat.ranges <- GenomicRanges::makeGRangesFromDataFrame(pat, 
                                                        keep.extra.columns = T, 
                                                        ignore.strand = T, 
                                                        end.field = 'end', 
                                                        start.field = 'start')
  
  pat.filt <- filter_pat(pat.ranges)

  # Find overlaps, taking into account overlap
  # Find all the overlaps, including anything with overhangs.
  # minoverlap by default is set to 3.
  marker.pat.overlaps.any <- findOverlaps(query = pat.filt, 
                                      subject = marker.ranges, 
                                      minoverlap = 3, type = "any")
  
  if(length(marker.pat.overlaps.any)==0){
    stop("ERROR: PAT does not have any overlaps with marker!")
  }
  
  # Find overlaps of "ideal scenarios", where the entire pat read is within the marker region
  marker.pat.overlaps.within <- findOverlaps(query = pat.filt, 
                                             subject = marker.ranges, 
                                             minoverlap = 3, type = "within")
  
  if(isFALSE(split_reads)){
    warning("split_read = FALSE; Returning only overlaps of PAT reads within marker regions.")
    return(list(overlaps = marker.pat.overlaps.within, pat.gr = pat.filt, marker.gr = marker.ranges))
  }
  
  # Remove the ideal scenarios, leaving only those with overhangs
  # Trim and split the overhanging reads in the next loop
  marker.pat.overlaps.notwithin <- marker.pat.overlaps.any[-S4Vectors::match(x = marker.pat.overlaps.within, table = marker.pat.overlaps.any)]
  
  # Split aligned, overhanging reads from pat file by CpG index
  # First, get indices of all overhanging reads
  pat_ind_all <- marker.pat.overlaps.notwithin@from
  
  # Then, loop through overhanging read indices
  # pat.overhanging.dt <- as.data.table(pat)[pat_ind_all,]
  # markers.overhanging <- marker.ranges[pat_ind_all,]
  # substr(pat.dt$read)
  # 
  
  
  
  
  res <- pbapply::pblapply(pat_ind_all, function(pat_ind){
    # Get pat grange object that corresponds to overhanging pat read
    pat.granges_i <- pat.ranges[pat_ind]
    # Find marker index that corresponds with overhanging pat read
    marker_ind <- marker.pat.overlaps.any@to[marker.pat.overlaps.any@from == pat_ind]
    # Get the marker granges object that corresponds with marker_ind
    marker.granges_i <- marker.ranges[marker_ind]
    pat.ranges.to_add <- GenomicRanges::GRanges()
    #pat.ranges.to_add <- list()
    new_ind <- 1
    for(j in 1:length(marker.granges_i)){
      marker.granges_ij <- marker.granges_i[j]
      new_grange <- GenomicRanges::intersect(pat.granges_i, marker.granges_ij)
      rel_start <- new_grange@ranges@start - pat.granges_i@ranges@start
      rel_end <- rel_start + new_grange@ranges@width
      new_grange$read <- substring(pat.granges_i$read, first = rel_start, last = rel_end)
      new_grange$nobs <- pat.granges_i$nobs
      pat.ranges.to_add[new_ind] <- new_grange
      # actual_start <- max(marker.ranges_ij@start, pat.ranges_i@start)
      # rel_start <- (actual_start - pat.ranges_i@start) + 1
      # actual_end <- min(marker.ranges_ij@start+(marker.ranges_ij@width-1), pat.ranges_i@start+(pat.ranges_i@width-1))
      # rel_end <- rel_start + (actual_end - actual_start) + 1
      # pat.ranges.to_add[new_ind] <- GenomicRanges::GRanges(
      #   seqnames = pat.granges_i@seqnames,
      #   ranges = IRanges(start = actual_start, end = actual_end),
      #   read = substring(pat.granges_i$read, first = rel_start, last = rel_end),
      #   nobs = pat.granges_i$nobs
      # )
      new_ind = new_ind + 1
    }
    return(pat.ranges.to_add)
  }, cl = n_threads)
  
  # Generate new pat, filter, and re-align
  res_lengths <- sapply(res, length)
  if(any(res_lengths==0)){
    res.grl <- GRangesList(res[res_lengths > 0])
    pat.ranges.split_reads <- unlist(res.grl)
  }else{
    pat.ranges.split_reads <- do.call(c, res)
  }
  #This still returns a list of granges objects...
  pat.ranges.all <- unlist(c(pat.ranges.split_reads, pat.ranges[marker.pat.overlaps.within@from]))
  
  pat.filt <- filter_pat(pat = pat.ranges.all)
  
  # Second overlap
  marker.pat.overlaps.2 <- findOverlaps(query = pat.filt, 
                                      subject = marker.ranges, 
                                      minoverlap = 1, type = "any")
  
  return(list(overlaps = marker.pat.overlaps.2, pat.gr = pat.filt, marker.gr = marker.ranges))
}
