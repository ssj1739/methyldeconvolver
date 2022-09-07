
#' read_pat
#'
#' @param path character indicating path to PAT file. May be gzipped.
#' @param linelimit numeric. Default Inf.
#' @param verbose logical. Default FALSE.
#'
#' @return data.frame containing contents of pat file.
#' @export
#'
#' @examples
#' \dontrun{
#' read_pat(path = "path/to/pat_file.pat.gz")
#' }
read_pat <- function(path="data/Hep_all.pat.gz", linelimit = Inf, verbose = F){
  require(data.table)
  require(dplyr)
  if(verbose) message("Starting to read file.")
  pat = data.table::fread(path, nrows = linelimit, header = F)
  colnames(pat)[1:4] <- c("chr", "start", "read", "nobs")
  # Filter for read requirements:
  # - Must have at least 1 C or T
  # - Must be of length 3 or more
  if(verbose)
    message("Finished reading, now filtering.")
  pat.filt <- pat %>%
    dplyr::filter(nchar(read) >= 3) %>%
    dplyr::filter(grepl("C|T", read))
  
  if(verbose){
    message("Finished filtering.")
    message(paste0("Filtered out ", nrow(pat) - nrow(pat.filt), " rows"))
  }
    
  
  return(pat.filt)
}


#' read_pat2
#' Meant to be a more efficient way to stream in large pat files. Still not ready for export.
#' @param path 
#' @param chunksize 
#'
#' @return
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
#' @param path 
#' @param linelimit 
#'
#' @return data.frame containing marker file contents
#' @export
#'
#' @examples
#' \dontrun{
#' read_marker(path = "inst_data/human_mixintest_top25.txt")
#' }
read_marker <- function(path="data/Human_mixintest_top25.txt", linelimit = Inf){
  require(data.table)
  marker = data.table::fread(path, nrows = linelimit, header = T)
  
  # TODO: Valudate marker file format
  
  return(marker)
}
