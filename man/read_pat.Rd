% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_files.R
\name{read_pat}
\alias{read_pat}
\title{read_pat}
\usage{
read_pat(path = "data/ref/Hep_all.pat.gz", linelimit = Inf, verbose = F, ...)
}
\arguments{
\item{path}{character indicating path to PAT file. May be gzipped.}

\item{linelimit}{numeric. Default Inf.}

\item{verbose}{logical. Default FALSE.}
}
\value{
data.frame containing contents of pat file.
}
\description{
Reads in PAT-format files (output from wgbstools).
Also includes
}
\examples{
\dontrun{
read_pat(path = "path/to/pat_file.pat.gz")
}
}
