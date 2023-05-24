# Analysis of in-silico mix-in data
library(methyldeconvolveR)

reference <- learn_reference(marker.file = "~/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/markers/markers.ALL.bulk-immune.top200.q0.05.margin0.3.bed",
                            pat.dir = "~/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/pat/",
                            verbose = T, n_threads = 6)
reference <- readRDS("~/Tools/methyldeconvolver/data/reference_4-10-23.rds")

sample_pat <- "~/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/test_mixin/healthy.test_1.pat.gz"

results <- deconvolute_sample(sample_pat = sample_pat, reference = reference, quiet = F, n_threads = 6)
