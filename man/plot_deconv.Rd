% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_deconv.R
\name{plot_deconv}
\alias{plot_deconv}
\title{plot_deconv}
\usage{
plot_deconv(results, colormap, ...)
}
\arguments{
\item{results}{\itemize{
\item either the simple output of deconvolute_sample, or the
"best_result" slot in the full output of deconovlute_sample
}}

\item{x.lab}{Character passed to ggplot2::xlab()}

\item{y.lab}{Character passed to ggplot2::ylab()}

\item{...}{Other arguments passed on to ggplot2::theme().}
}
\description{
Quick visualization of deconvolution results from a single sample
}
\examples{
\dontrun{
  results <- deconvolute_sample("path/to/sample.pat.gz")
  plot_deconv(results)
}
}
