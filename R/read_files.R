
#' read_pat
#'
#' @param path character indicating path to PAT file. May be gzipped.
#' @param linelimit numeric. Default Inf.
#' @param verbose logical. Default FALSE.
#' 
#' @description
#' Reads in PAT-format files (output from wgbstools). 
#' Also includes 
#'
#' @return data.frame containing contents of pat file.
#' @export
#'
#' @examples
#' \dontrun{
#' read_pat(path = "path/to/pat_file.pat.gz")
#' }
read_pat <- function(path="data/ref/Hep_all.pat.gz", 
                     linelimit = Inf,
                     verbose = F,
                     ...){
  require(data.table)
  require(dplyr)
  if(verbose) message("Starting to read file.")
  pat = data.table::fread(file = path, nrows = linelimit, header = F)
  colnames(pat)[1:4] <- c("chr", "start", "read", "nobs")
  
  # Filter for read requirements:
  # - Must have at least 1 C or T
  # - Must be of length 3 or more
  if(verbose)
    message("Finished reading, now filtering.")
  
  pat.filt <- filter_pat(pat, ...)
  
  return(pat.filt)
}


#' read_pat2
#' Meant to be a more efficient way to stream in large pat files. Still not ready for export.
#' @param path 
#' @param chunksize 
#'
#' @return data.frame containing contents of pat file
#'
#' @examples
#' \dontrun{
#' read_pat2(path = "path/to/pat_file.pat.gz")
#' }
read_pat2 <- function(path="data/Hep_all.pat.gz", chunksize = 5000){
  require(utils)
  # First check number of lines
  if(grepl("gz", path)){
    con <- gzfile(path, open = "r")
  }else{
    con <- file(description = path, open = "r")
  }
  datachunk <- utils::read.table(con, nrows = chunksize)
  # Stream in file, chunk by chunk - figure out some way to process each chunk here?
  
  on.exit(close(con))
}


#' read_marker
#'
#' @param path character indicating path to PAT file. May be gzipped.
#' @param linelimit numeric. Default Inf.
#'
#' @return data.frame containing marker file contents
#' @export
#'
#' @examples
#' \dontrun{
#' read_marker(path = "inst_data/human_mixintest_top25.txt")
#' }
read_marker <- function(path="data/Human_mixintest_top25.txt", no_reduction = F){
  require(data.table)
  require(GenomicRanges)
  require(tidyverse)
  
  marker = data.table::fread(path, nrows = Inf, header = F)
  if(!grepl("[1-9]",marker[1,1])){
    marker <- marker[2:nrow(marker),]
  }
  colnames(marker)[1:8] <- c("chr", "startChr", "endChr", "startCpG", "endCpG", "target", "name", "direction")

  # Normalize case of target cell types
  marker$target <- tolower(marker$target)
  
  # Merge overlapping marker regions
  marker.gr <- GenomicRanges::makeGRangesFromDataFrame(marker, 
                                                    start.field = "startCpG", 
                                                    end.field = "endCpG",
                                                    keep.extra.columns = T)
  
  if(isTRUE(no_reduction)){
    return(marker.gr)
  }
    
  
  marker.gr.red <- GenomicRanges::reduce(marker.gr)

  reannotate.ol <- GenomicRanges::findOverlaps(marker.gr.red, marker.gr)
  reannotate.target <- character(length(marker.gr.red))
  reannotate.oi <- character(length(marker.gr.red))
  reannotate.n_merged <- numeric(length(marker.gr.red))
  for(i in 1:length(marker.gr.red)){
    original.ind <- reannotate.ol@to[reannotate.ol@from==i]
    reannotate.target[i] <- paste0(unique(marker.gr$target[original.ind]), collapse = ",")
    reannotate.oi[i] <- paste0(original.ind, collapse = ',')
    reannotate.n_merged[i] <- length(original.ind)
  }
  
  marker.gr.red$target <- reannotate.target
  marker.gr.red$original.ind <- reannotate.oi
  marker.gr.red$n_merged <- reannotate.n_merged
  
  
  out_marker <- marker.gr.red
  return(out_marker)
}

