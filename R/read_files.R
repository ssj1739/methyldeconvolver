
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
#' @importFrom data.table fread
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
  if(verbose) message("Starting to read file.")
  pat = data.table::fread(file = path, nrows = linelimit, header = F)
  colnames(pat)[1:4] <- c("chr", "start", "read", "nobs")
  
  # Filter for read requirements:
  # - Must have at least 1 C or T
  # - Must be of length 3 or more
  if(verbose)
    message("Finished reading, now filtering.")
  
  pat.filt <- filter_pat(pat, verbose = verbose, ...)
  
  return(pat.filt)
}


#' read_pat2
#' Meant to be a more efficient way to stream in large pat files. Still not ready for export.
#' @param path 
#' @param chunksize 
#' 
#' @importFrom utils read.table
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
  
  return(datachunk)
  on.exit(close(con))
}


#' read_marker
#'
#' @param path character indicating path to PAT file. May be gzipped.
#' @param linelimit numeric. Default Inf.
#'
#' @return data.frame containing marker file contents
#' @importFrom utils read.table 
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#'
#' @examples
#' \dontrun{
#' # Access path to marker using system.file(package = "methyldeconvolveR", "extdata", "marker_example_top100.txt")
#' read_marker(path = "path/to/marker.txt")
#' }
read_marker <- function(path="", no_reduction = F, header = T, sep = "\t"){
  require(utils)
  require(GenomicRanges)

  marker = utils::read.table(file = path, header = header, sep = sep, comment.char = "#")
  if(!grepl("[1-9]",marker[1,1])){
    marker <- marker[2:nrow(marker),]
  }
  if(!"label" %in% colnames(marker)){
    # Assume default marker output
    # TODO: Assess if these colnames would be correct somehow?
    colnames(marker)[1:6] <- c("chr", "startChr", "endChr", "startCpG", "endCpG", "label")
  }
  if(any(c("start", "end") %in% colnames(marker))){
    which.start <- which(colnames(marker) %in% c("start"))
    which.end <- which(colnames(marker) %in% c("end"))
    colnames(marker)[which.start] <- "startChr"
    colnames(marker)[which.end] <- "endChr"
  }

  # Normalize case and punctuation of target cell types
  marker$label <- marker$label %>% tolower() %>% 
    gsub(pattern = "[[:punct:]]", replacement = "") %>%
    gsub(pattern = " ", replacement = "")
  
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
  reannotate.p.value <- reannotate.q.value <- numeric(length(marker.gr.red))
  
  pb <- pbapply::startpb(min = 0, max = length(marker.gr.red))
  for(i in 1:length(marker.gr.red)){
    original.ind <- reannotate.ol@to[reannotate.ol@from==i]
    reannotate.target[i] <- paste0(unique(marker.gr$label[original.ind]), collapse = ",")
    reannotate.oi[i] <- paste0(original.ind, collapse = ',')
    reannotate.n_merged[i] <- length(original.ind)
    if("p.value" %in% colnames(mcols(marker.gr)))
      reannotate.p.value[i] <- paste0(marker.gr$p.value[original.ind], collapse = ",")
    if("q.value" %in% colnames(mcols(marker.gr)))
      reannotate.q.value[i] <- paste0(marker.gr$q.value[original.ind], collapse = ",")
    pbapply::setpb(pb, value = i)
  }
  
  marker.gr.red$label <- reannotate.target
  marker.gr.red$original.ind <- reannotate.oi
  marker.gr.red$n_merged <- reannotate.n_merged
  
  if("p.value" %in% colnames(mcols(marker.gr)))
    marker.gr.red$p.val <- reannotate.p.value
  if("q.value" %in% colnames(mcols(marker.gr)))
    marker.gr.red$q.val <- reannotate.q.value
  
  return(marker.gr.red)
}

