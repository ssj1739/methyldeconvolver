#' use_reference
#'
#' @description Function to access the published reference atlas described in McDeed et al.
#' 
#' @return reference object
#' @export
#' @examples
#' 
#' ref <- use_reference()
#' names(ref)
#' 
#' 
use_reference <- function(){
  data("reference", package = "methyldeconvolveR", envir = environment())
  return(reference)
}
