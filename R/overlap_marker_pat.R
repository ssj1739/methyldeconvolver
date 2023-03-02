
#' Overlap marker and APT files
#' 
#'@description Merging the PAT and Marker files together
#'@param pat Some description of PAT files
#'@param marker Some description of Marker files
#'@param n_threads Number of threads to use (default = 1)
#'
#'@export
#'
overlap_marker_pat <- function(pat, marker, n_threads = 1){
  require(GenomicRanges)
  require(dplyr)
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

  
  # Find overlaps, taking into account overlap
  marker.pat.overlaps <- findOverlaps(query = pat.ranges, 
                                      subject = marker.ranges, 
                                      minoverlap = 1, type = "any")
  
  # Split aligned reads from pat file by CpG index
  # multiple_marker_alignment_i <- table(marker.pat.overlaps@from)
  # multiple_marker_alignment <- as.numeric(names(multiple_marker_alignment_i)[which(multiple_marker_alignment_i > 1)])
  # pat.ranges.removed_multiple_aligned <- pat.ranges[-multiple_marker_alignment]
  # pat.ranges.to_add <- GenomicRanges::GRanges()

  pat_ind_all <- marker.pat.overlaps@from
    
  res <- pbapply::pbsapply(pat_ind_all, function(pat_ind){
    pat.granges_i <- pat.ranges[pat_ind]
    pat.ranges_i <- pat.ranges[pat_ind]@ranges
    marker_ind <- marker.pat.overlaps@to[marker.pat.overlaps@from == pat_ind]
    marker.ranges_i <- marker.ranges[marker_ind]@ranges
    pat.ranges.to_add <- GenomicRanges::GRanges()
    new_ind <- 1
    for(j in 1:length(marker.ranges_i)){
      marker.ranges_ij <- marker.ranges_i[j]
      actual_start <- max(marker.ranges_ij@start, pat.ranges_i@start)
      rel_start <- (pat.ranges_i@start - actual_start) + 1
      actual_end <- min(marker.ranges_ij@start+(marker.ranges_ij@width-1), pat.ranges_i@start+(pat.ranges_i@width-1))
      rel_end <- actual_end - actual_start + 1
      pat.ranges.to_add[new_ind] <- GenomicRanges::GRanges(
        seqnames = pat.ranges@seqnames[pat_ind], 
        ranges = IRanges(start = actual_start, end = actual_end),
        read = substring(pat.granges_i$read, first = rel_start, last = rel_end),
        nobs = pat.granges_i$nobs
      )
      new_ind = new_ind + 1
    }
    return(pat.ranges.to_add)
  }, cl = n_threads)
  # Generate new pat file and re-align
  pat.ranges.split_reads <- do.call(`c`, res)
  
  # Second overlap
  marker.pat.overlaps.2 <- findOverlaps(query = pat.ranges.split_reads, 
                                      subject = marker.ranges, 
                                      minoverlap = 1, type = "any")
  
  
  return(list(overlaps = marker.pat.overlaps.2, pat.gr = pat.ranges.split_reads, marker.gr = marker.ranges))
}
