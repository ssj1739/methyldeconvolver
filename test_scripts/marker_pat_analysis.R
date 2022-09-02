library(methyldeconvolveR)
library(GenomicRanges)
library(stringr)
library(magrittr)
library(tidyverse)


marker <- read_marker("data/Human_mixintest_top25.txt")
pat <- read_pat("data/Test25_1.pat.gz")
#colnames(pat) <- c("chr", "start", "read", "nobs")
#pat$end <- pat$start + nchar(pat$read)

ref <- learn_reference(marker.file = "data/Human_mixintest_top25.txt", pat.dir = "data/ref/", save.output = "test_reference.rds", verbose = T)

res_unweighted <- deconvolute_sample(sample.pat.path = "data/Test25_1.pat.gz", reference = ref, verbose = T, n_threads = 6, num_of_inits = 10)


res_weighted <- deconvolute_sample_weighted(sample.pat.path = "data/Test25_1.pat.gz", reference = ref, verbose = T, n_threads = 6, num_of_inits = 10)

# Calculating the number of read * observations for each marker
aligned.read.nobs <- overlaps.list$pat.gr$nobs[overlaps.list$overlaps@from]
table(rep(overlaps.list$overlaps@to, aligned.read.nobs))


# marker.ranges <- GenomicRanges::makeGRangesFromDataFrame(marker[,-c("start", "end")], start.field = "startCpG", end.field = "endCpG", keep.extra.columns = T)
# pat.ranges <- GenomicRanges::makeGRangesFromDataFrame(pat, keep.extra.columns = T)
# 
# marker.pat.overlaps <- findOverlaps(query = pat.ranges, subject = marker.ranges, minoverlap = 1, type = "any")
# 
# beta_fits <- sapply(unique(marker.pat.overlaps@to), function(x){
#   pat.index <- marker.pat.overlaps@from[marker.pat.overlaps@to == x]
#   if(length(pat.index)==0){
#     return(NULL)
#   }
#   # Potentially add filter step to remove <=2 length
# 
#   pat.subset <- pat.ranges[pat.index]
#   nTs = stringr::str_count(pat.subset$read, "T")
#   nCs = stringr::str_count(pat.subset$read, "C")
#   meth.fraction <-  nCs / (nCs + nTs)
#   rep.meth.fraction <- rep(meth.fraction, pat.subset$nobs)
# 
#   # If there are not at least 3 distinct reads that map to a given site,
#   # and not at least 2 unique values (not all 0s or all 1s), then consider modeling the marker:
#   
#   fit <- NA
#   if(length(rep.meth.fraction) >= 3 & length(unique(rep.meth.fraction)) >= 2){
#     # First try mle:
#     fit = try(fitdistrplus::fitdist(data = rep.meth.fraction, distr = "beta", method = "mle"), silent = T)
#     
#     # If mle doesn't work, then try mme:
#     if(class(fit)=="try-error"){ 
#       fit = fitdistrplus::fitdist(data = rep.meth.fraction, distr = "beta", method = "mme")
#     }
#   }
#   return(fit)
# }) %>% magrittr::set_names(unique(marker.pat.overlaps@to))
# 
# 
# 
# 
# beta_fits <- list()
# for(x in unique(marker.pat.overlaps@to)){
#   pat.index <- marker.pat.overlaps@from[marker.pat.overlaps@to == x]
#   if(length(pat.index)==0){
#     return(NULL)
#   }
#   # Potentially add filter step to remove <=2 length
#   
#   pat.subset <- pat.ranges[pat.index]
#   nTs = stringr::str_count(pat.subset$read, "T")
#   nCs = stringr::str_count(pat.subset$read, "C")
#   meth.fraction <-  nCs / (nCs + nTs)
#   rep.meth.fraction <- rep(meth.fraction, pat.subset$nobs)
#   
#   # If there are not at least 3 distinct reads that map to a given site,
#   # and not at least 2 unique values (not all 0s or all 1s), then consider modeling the marker:
#   if(length(rep.meth.fraction) < 3 | length(unique(rep.meth.fraction)) < 2){
#     beta_fits[[x]] <- NA
#   }else{
#     # First try mle:
#     fit = try(fitdistrplus::fitdist(data = rep.meth.fraction, distr = "beta", method = "mle"), silent = T)
#     
#     # If mle doesn't work, then try mme:
#     if(class(fit)=="try-error"){ 
#       fit = fitdistrplus::fitdist(data = rep.meth.fraction, distr = "beta", method = "mme")
#     }
#     # Save fit output:
#     beta_fits[[x]] <- fit$estimate
#   }
# }

# Fit custom beta dist
beta_ll <- function(theta, x){
  N <- length(x)
  alpha <- theta[1]
  beta <- theta[2]
  beta_ll <- ((alpha-1)*sum((log(x)))) + ((beta-1)*sum(log(1-x))) - N*log(beta(alpha, beta))
  return(-beta_ll)
}

test.vec <- c(rep(1, 10), rbeta(n = 10, shape = 0.1, shape2 = 10))

res = optim(par = c(0.01, 0.01), fn = beta_ll, x = test.vec, method = "L-BFGS-B", lower = 1e-3, upper = 100)

res$par

# Computing width of overlap:
overlaps <- pintersect(overlaps.list$pat.gr[queryHits(overlaps.list$overlaps)], overlaps.list$marker.gr[subjectHits(overlaps.list$overlaps)])
percentOverlap <- width(overlaps) / width(overlaps.list$marker.gr[subjectHits(overlaps.list$overlaps)])
overlaps.list$overlaps.over50 <- overlaps.list$overlaps[percentOverlap > 0.5]

# Which reads overlap with more than 1 marker?
multiple_overlaps.pat.gr <- overlaps.list$pat.gr[as.numeric(names(which(table(overlaps.list$overlaps@from)==2)))]
# Realigned those reads that had multiple marker overlaps (because subsetting changed indices):
hits <- findOverlaps(multiple_overlaps.pat.gr, overlaps.list$marker.gr)
overlaps <- pintersect(multiple_overlaps.pat.gr[queryHits(hits)], overlaps.list$marker.gr[subjectHits(hits)])
multiple_overlaps.pat.gr
