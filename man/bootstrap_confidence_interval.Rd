% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_confidence_interval.R
\name{bootstrap_confidence_interval}
\alias{bootstrap_confidence_interval}
\title{bootstrap_confidence_interval}
\usage{
bootstrap_confidence_interval(
  num_boots = 1000,
  sample_pat,
  reference,
  n_threads = 1
)
}
\arguments{
\item{num_boots}{number of bootstraps}

\item{sample_pat}{pat file to bootstrap from}

\item{reference}{reference object}

\item{n_threads}{number of threads (default 1)}
}
\value{
confidence intervals around each proportion estimate
}
\description{
bootstrap_confidence_interval
}
\examples{
\dontrun{
bootstrap_confidence_interval(num_boots = 1000, sample_pat = "sample.pat.gz", reference = ref)
}

}
