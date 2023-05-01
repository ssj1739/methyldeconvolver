# Analysis of in-silico mix-in data
library(methyldeconvolveR)

reference <- readRDS("~/Tools/methyldeconvolver/data/reference_4-10-23.rds")

sample_pat <- "~/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/test_mixin/healthy.test_1.pat.gz"

results <- deconvolute_sample(sample_pat = sample_pat, reference = reference, quiet = F, n_threads = 6)
