# Full test
library(methyldeconvolveR)
library(tidyverse)

#ref <- learn_reference(marker.file = "data/FullTest/Sid_testmarkers.bed", pat.dir = "data/FullTest/ref/", save.output = "data/FullTest/FullTest_ref.rds")
ref <- readRDS("data/FullTest/FullTest_ref.rds")
sample.pat <- read_pat(path = "data/FullTest/test_mix_1.pat.gz", verbose = T)

true_distr <- rep(sample.pat$V5, times = sample.pat$nobs)
cell_types <- unique(sample.pat$V5)

count_distr <- numeric(5)
names(count_distr) <- cell_types
for(i in 1:length(cell_types)){
  print(paste0("Completing ", cell_types[i]))
  count_distr[i] <- sum(true_distr == cell_types[i])
}

true_prop <- count_distr / sum(count_distr)
true_prop <- true_prop[sort(names(true_prop))]
print(true_prop)

res_unweighted <- deconvolute_sample(sample.pat.path = "data/FullTest/test_mix_1.pat.gz", retain_alphas = T, reference = ref, verbose = T, n_threads = 6, num_of_inits = 100)
res_weighted <- deconvolute_sample_weighted(sample.pat.path = "data/FullTest/test_mix_1.pat.gz", reference = ref, verbose = T, n_threads = 6, num_of_inits = 100)


lapply(res_unweighted, function(x){
  caret::RMSE(x$last_alpha, true_prop[names(x$last_alpha)])
})

weighted_iter_df <- lapply(res_weighted, function(x){
  return(x$last_alpha)
}) %>%
  bind_rows() %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "cell_type", values_to = "estimated_proportion") %>%
  mutate(Estimate = "weighted")

unweighted_iter_df <- lapply(res_unweighted, function(x){
  return(x$last_alpha)
}) %>%
  bind_rows() %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "cell_type", values_to = "estimated_proportion") %>%
  mutate(Estimate = "unweighted")

truth_df <- data.frame(
  cell_type = names(true_prop),
  estimated_proportion = true_prop,
  Estimate = "truth"
)

combined_iter_df <- bind_rows(weighted_iter_df, unweighted_iter_df, truth_df)

ggplot(combined_iter_df, aes(x = cell_type, y = estimated_proportion, color = Estimate)) +
  geom_boxplot()

ggsave(filename = "figs/MixIn_WeightedvsUnweighted_Boxplot.png", width = 8, height = 5, dpi = 300)

lapply(res_weighted, function(x){
  caret::RMSE(x$last_alpha, true_prop[names(x$last_alpha)])
})

# Visualize results
res_df <- data.frame(
  cell_types = names(true_prop),
  truth = true_prop,
  unweighted = res_unweighted[[1]]$last_alpha,
  weighted = res_weighted[[1]]$last_alpha
)

res_df_clean <- res_df %>%
  tidyr::pivot_longer(cols = -cell_types, names_to = "Estimate", values_to = "Proportions")

ggplot(res_df_clean, aes(x = cell_types, y = Proportions, color = Estimate)) +
  geom_boxplot()

ggsave(filename = "figs/MixIn_WeightedvsUnweighted.png", width = 8, height = 5, dpi = 300)

mad_unweighted_df = as.data.frame(res_unweighted[[1]]$iter_mad) %>%
  magrittr::set_colnames("MAD_value") %>%
  mutate(iteration = 1:nrow(.))

ggplot(mad_unweighted_df, aes(x = iteration, y = MAD_value)) +
  geom_point()

ggsave("mad_unweighted.png", width = 7, height = 4)


mad_weighted_df = as.data.frame(res_weighted[[1]]$iter_mad) %>%
  magrittr::set_colnames(c("mad.1", "mad.0")) %>%
  mutate(iteration = 1:nrow(.)) %>%
  pivot_longer(cols = c(mad.1, mad.0), names_to = "MAD_label", values_to = "MAD_value")

ggplot(mad_weighted_df, aes(x = iteration, y = MAD_value, color = MAD_label)) +
  geom_point() +
  facet_wrap(~MAD_label, nrow = 2)

ggsave("mad_weighted.png", width = 7, height = 7)
