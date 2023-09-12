#' plot_deconv
#'
#' @usage plot_deconv(results, colormap, ...)
#' @description
#' Quick visualization of deconvolution results from a single sample
#' 
#' 
#' @param results - either the simple output of deconvolute_sample, or the 
#' "best_result" slot in the full output of deconovlute_sample
#' @param x.lab Character passed to ggplot2::xlab()
#' @param y.lab Character passed to ggplot2::ylab() 
#' @param ... Other arguments passed on to ggplot2::theme().
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @export
#'
#' @examples
#' \dontrun{
#'   results <- deconvolute_sample("path/to/sample.pat.gz")
#'   plot_deconv(results)
#' }
plot_deconv <- function(results, x.lab = "Sample", y.lab = "Proportions", ...){
  require(ggplot2)
  
  plot_df <- data.frame(
    sample_name = "Sample",
    cell_types = names(results),
    proportions = results
  )
  
  colours_custom = c("#00A6FB","#F8A079","#0582CA","#F4743B","#70AE6E","#006494", "#E79345","#483C46","#D9B14F", "#42555C","#CCD059", "#3C6E71","#BEEE62","#568E70")
  
  plot <- ggplot(data = plot_df, aes(x = sample_name, y = proportions, fill = cell_types)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colours_custom) +
    theme_minimal(...) +
    xlab(x.lab) +
    ylab(y.lab)
  
  return(plot)
}
