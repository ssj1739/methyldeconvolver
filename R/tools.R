#' encode_binary
#'
#' @param read 
#'
#' @return
#' @export
#'
#' @examples
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
