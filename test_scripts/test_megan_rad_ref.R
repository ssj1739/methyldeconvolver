# Learn Full Reference
library(methyldeconvolveR)

#pat.dir <- "~/Documents/Research/Wellstein/Projects/Megan_Radiation/Methylomes/PAT/"
pat.dir <- "~/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/pat/"
#marker.file <- "~/Documents/Research/Wellstein/Projects/Megan_Radiation/Methylomes/Marker/markers.ALL.bed"
marker.file <- "~/Data/TestMixIn/from_megan_rad/Markers/FinalMarkers/markers.ALL.bed"

marker <- read_marker(marker.file)
pat.files = dir(pattern = "*.pat.gz$", pat.dir, full.names = F)

reference <- learn_reference(marker.file = marker.file, pat.dir = pat.dir, verbose = T, n_threads = 12, split_reads = F)
saveRDS(reference, file = "~/Tools/methyldeconvolver/data/reference_4-21-23.rds")
# reference <- readRDS("~/Tools/methyldeconvolver/data/reference_4-10-23.rds")

# 
# pat.files <- dir(pat.dir,pattern = "*.pat.gz$", full.names = T)
# # 
#  result_deconv_self <- list()
#  pb <- pbapply::startpb(min = 0, max = length(pat.files))
#  for(i in seq_along(pat.files)){
#    pf <- pat.files[i]
#    res <- deconvolute_sample(sample_pat = pf, reference = reference, n_threads = 4, num_of_inits = 1)
#    result_deconv_self[[pf]] <- res$last_alpha
#    pbapply::setpb(pb, value = i)
#  }
#  closepb(pb)
# # saveRDS(result_deconv_self, file = "All_results_full_ref.rds")
# 
# result_deconv_self <- readRDS("data/megan_rad/All_results_full_ref.rds")
# 
# pat.cell_types <- sapply(pat.files, function(x){
#   s = strsplit(x, split = "/", fixed = T)[[1]]
#   return(strsplit(s[length(s)], split = "_", fixed = T)[[1]][1])
# })
# result_df <- bind_rows(result_deconv_self, .id = "PAT file")
# result_df$target_cell_type <- pat.cell_types[1:nrow(result_df)]
# 
# g1 <- heatmaply::ggheatmap(result_df %>% select(-`PAT file`, -target_cell_type),
#                            row_side_colors = tolower(result_df$target_cell_type),
#                      Rowv = F, Colv = F, 
#                      dendrogram = "row",
#                      showticklabels = c(T, F))
# ggsave(plot = g1, "figs/SelfDeconvolution_heatmap_raw.png", width = 12, height = 6, dpi = 300)
# 
# g2 <- heatmaply::ggheatmap(result_df %>% select(-`PAT file`),
#                            showticklabels = c(T, F))
# 
# ggsave(plot = g2, "figs/SelfDeconvolution_heatmap_clust.png", width = 12, height = 6, dpi = 300)
# 
# markers_per_celltype <- as.data.frame(table(reference$marker$target))
# 
# summarized_result_df <- result_df %>%
#   tidyr::pivot_longer(cols = -c("PAT file", "target_cell_type"), names_to = "deconv_cell_type", values_to = "proportion_estimate") %>%
#   group_by(deconv_cell_type, target_cell_type) %>%
#   summarize(mean(proportion_estimate)) %>%
#   filter(tolower(target_cell_type)==deconv_cell_type) %>%
#   merge(., markers_per_celltype, by.x = "deconv_cell_type", by.y = "Var1")
# 
# ggplot(data = summarized_result_df, aes(y = `mean(proportion_estimate)`, x = Freq, color = target_cell_type)) +
#   geom_point() +
#   labs(y = "Avg. Proportion Estimate for Self-Deconvolution", x = "Number of Markers for Cell-type" ) 
# 
# ggsave(filename = "figs/SelfDeconvolution_Markers_vs_Proportions.png")
# 
# 
# heatmaply::heatmaply(result_df %>% select(-`PAT file`))
