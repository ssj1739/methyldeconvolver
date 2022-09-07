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
