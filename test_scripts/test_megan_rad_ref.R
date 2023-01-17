# Learn Full Reference
library(methyldeconvolveR)

pat.dir <- "~/Documents/Research/Wellstein/Projects/Megan_Radiation/Methylomes/PAT/"
marker.file <- "~/Documents/Research/Wellstein/Projects/Megan_Radiation/Methylomes/Marker/markers.ALL.bed"

#reference <- learn_reference(marker.file = marker.file, pat.dir = pat.dir, verbose = T)
reference <- readRDS("complete_reference_megan_rad_1-11-23.rds")

pat.files <- dir(pat.dir, full.names = T)

# result_deconv_self <- list()
# for(pf in pat.files){
#   res <- deconvolute_sample(sample_pat = pf, reference = reference, n_threads = 1, num_of_inits = 1)
#   result_deconv_self[[pf]] <- res$last_alpha
# }
# saveRDS(result_deconv_self, file = "All_results_full_ref.rds")

result_deconv_self <- readRDS("All_results_full_ref.rds")

pat.cell_types <- sapply(pat.files, function(x){
  s = strsplit(x, split = "/", fixed = T)[[1]]
  return(strsplit(s[length(s)], split = "_", fixed = T)[[1]][1])
})
result_df <- bind_rows(result_deconv_self, .id = "PAT file")
result_df$target_cell_type <- pat.cell_types[1:nrow(result_df)]

markers_per_celltype <- as.data.frame(table(reference$marker$target))

summarized_result_df <- result_df %>%
  tidyr::pivot_longer(cols = -c("PAT file", "target_cell_type"), names_to = "deconv_cell_type", values_to = "proportion_estimate") %>%
  group_by(deconv_cell_type, target_cell_type) %>%
  summarize(mean(proportion_estimate)) %>%
  filter(tolower(target_cell_type)==deconv_cell_type) %>%
  merge(., markers_per_celltype, by.x = "deconv_cell_type", by.y = "Var1")

ggplot(data = summarized_result_df, aes(y = `mean(proportion_estimate)`, x = Freq, color = target_cell_type)) +
  geom_point() +
  labs(y = "Avg. Proportion Estimate for Self-Deconvolution", x = "Number of Markers for Cell-type" ) 

ggsave(filename = "figs/SelfDeconvolution_Markers_vs_Proportions.png")


heatmaply::heatmaply(result_df %>% select(-`PAT file`))
