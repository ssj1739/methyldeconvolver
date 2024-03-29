% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_files.R
\name{read_marker}
\alias{read_marker}
\title{read_marker}
\usage{
read_marker(path = "", no_reduction = F, header = T, sep = "\\t")
}
\arguments{
\item{path}{character indicating path to PAT file. May be gzipped.}

\item{linelimit}{numeric. Default Inf.}
}
\value{
data.frame containing marker file contents
}
\description{
read_marker
}
\examples{
\dontrun{
# Access path to marker using system.file(package = "methyldeconvolveR", "extdata", "marker_example_top100.txt")
read_marker(path = "path/to/marker.txt")
}
}
