
<!-- README.md is generated from README.Rmd. Please edit that file -->

# methyldeconvolveR

<!-- badges: start -->
<!-- badges: end -->

This package exists to allow users to perform cell-type deconvolution
using whole-genome bisulfite sequencing.

## Installation

You can install the development version of methyldeconvolveR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ssj1739/methyldeconvolver")
```

## Learn Reference

Learn a reference set from given reference pat files of cell types of interest. Requires path to a marker file, path to PAT file directory, and output directory.

Note: Name formats of the reference PAT files should follow the following convention:
  - Files should end in .pat.gz (bgzipped PAT files outputted from wgbstools).
  - Files should have the name of the cell type (which perfectly matches the name in the marker file)
    followed by an underscore, followed by any other names.
  - e.g. Bcell_Sample_Name.pat.gz

``` r
library(methyldeconvolveR)
## basic example code
learn_reference(marker.file = "marker.txt", pat.dir = "data/ref/", save.output = "reference.rds")
```

## Deconvolute Sample
Deconvolution algorithm for a single liquid biopsy sample. Requires a path to sample in PAT format and reference .rds object.
Default "simple" output is a named vector of cell-type proportion estimates that maximize the log-likelihood function.
EM initializations and stopping criteria can be adjusted by the user from defaults.
Multiple threads/cores can be used to parallelize computations by specifying --n_threads

``` r
ref = "/path/to/reference/reference.rds"
deconvolute_sample_weighted(sample_pat = "sample_to_deconvolute.pat.gz", reference = ref)
```


