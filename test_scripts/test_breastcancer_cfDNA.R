# Analyze breast cancer cfmeDNA
library(methyldeconvolveR)
library(tidyverse)

reference <- readRDS("~/Tools/methyldeconvolver/data/reference_margin0.4.ALL_0723.rds")
reference <- readRDS("~/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/reference/reference_0723/em/reference_margin0.3.top100_0723.rds")

sample_pats <- dir(path = "~/Data/breastcancer_SRP167041", pattern = ".pat.gz$", recursive = T, full.names = T)
sample_annot <- read_csv("~/Data/breastcancer_SRP167041/SraRunTable (2).txt")

## Clean annotations ##
sample_annot$Stage <- "Healthy"
sample_annot$Stage[grepl("early", sample_annot$Isolate)] <- "Early"
sample_annot$Stage[grepl("late", sample_annot$Isolate)] <- "Late"

sample_annot$Subtype <- sapply(sample_annot$Isolate, function(x){
  y = strsplit(x, split = "\\", fixed = T)[[1]][3]
  substr(y, start = 3, stop = nchar(y))
})

result <- pblapply(sample_pats, FUN = deconvolute_sample, reference = reference, quiet = F, n_threads = 1, cl = length(sample_pats))
 
# result <- list()
# for(sample in sample_pats){
#   result[[basename(sample)]] <- deconvolute_sample(sample_pat = sample, reference = reference, quiet = F, n_threads = 1)
# }

result.df <- as.data.frame(sapply(result, function(x) x$best_result))
sample_names <- sapply(basename(sample_pats), function(x) strsplit(x, split = "_")[[1]][1])
colnames(result.df) <- sample_names
result.df <- result.df %>% rownames_to_column(var = "cell.type")

annot.result.df <- result.df %>% pivot_longer(cols = -cell.type, names_to = "Sample", values_to = "Prop") %>%
  merge(., sample_annot, by.x = "Sample", by.y = "Run")

res_fig <- ggplot(data = annot.result.df, aes(x = Sample, y = Prop, fill = cell.type)) +
  geom_bar(stat = "identity", color = "black") +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = Sample, y = -1, fill = Subtype)) +
  ggpubr::theme_pubr(x.text.angle = 90)
 
ggsave(res_fig, filename = "figs/BrCa_res_deconv_test.pdf", width = 8, height = 10)
