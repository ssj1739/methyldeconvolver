# Format PAT files for input to CelFie

require(data.table)
require(pbapply)
require(GenomicRanges)
require(tidyverse)

# Read in marker file
markerfile <- readr::read_delim("data/FullTest/Sid_testmarkers.bed", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

# Get ranges of interest from marker file
marker_df <- markerfile %>%
  select(chrom = chr, start = startCpG, end = endCpG) 

marker_ranges <- GenomicRanges::makeGRangesFromDataFrame(marker_df)

# Read in reference PAT files
### Point to PAT directory
patdir <- dir(path = "data/FullTest/ref", full.names = T)

#pb <- pbapply::startpb(min = 0, max = length(filename)*length(marker_ranges))
output_list <- pblapply(patdir, function(filename, marker_ranges){
  celltype = strsplit(strsplit(filename, split = "/", fixed = T)[[1]][4], split = "_", fixed = T)[[1]][1]
  
  # Read in pat files and convert to GRanges for subsetting
  pat <- methyldeconvolveR::read_pat(path = filename)
  pat_ranges <- pat %>%
    mutate(END = start + nchar(read)) %>% # add end column
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # Initialize empty columns
  meth <- numeric(length = length(marker_ranges))
  depth <- numeric(length = length(marker_ranges))
  
  # Subset based on marker region, compute depth and # of methylated reads
  # pb <- pbapply::startpb(min = 1, max = length(marker_ranges))
  # for(i in 1:length(marker_ranges)){
  #  setpb(pb = pb, value = i)
  #   range <- marker_ranges[i]
  #   ol <- IRanges::subsetByOverlaps(pat_ranges, range)
  #   depth[i] = sum(ol$nobs)
  #   meth[i] = sum(ol$nobs[grep("C", ol$read)])
  # }
  
  # System.time()
  pat_ranges_filt <- IRanges::subsetByOverlaps(pat_ranges, marker_ranges)
   output <- pblapply(seq_along(marker_ranges), function(i){
     ol <- IRanges::subsetByOverlaps(pat_ranges_filt, marker_ranges[i])
     depth = sum(ol$nobs)
     meth = sum(ol$nobs[grep("C", ol$read)])
     return(data.frame(meth = meth, depth = depth))
   }, cl = 4) %>% bind_rows()
  
  #output = data.frame(meth = meth, depth = depth)
  colnames(output) <- paste0(celltype, "_", colnames(output))
  return(output)
}, marker_ranges = marker_ranges)


ref_merged <- bind_cols(marker_df, output_list)

saveRDS(ref_merged, file = "test_scripts/merged_pat_celfie_ref.rds")

# Read sample pat
samplepatfile <- "data/FullTest/test_mix_1.pat.gz"
sample_pat <- methyldeconvolveR::read_pat(samplepatfile)
sample_ranges <- sample_pat %>%
  mutate(END = start + nchar(read)) %>% # add end column
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

sample_ranges_filt <- IRanges::subsetByOverlaps(sample_ranges, marker_ranges)

sample_output <- pblapply(seq_along(marker_ranges), function(i){
  ol <- IRanges::subsetByOverlaps(sample_ranges_filt, marker_ranges[i])
  depth = sum(ol$nobs)
  meth = sum(ol$nobs[grep("C", ol$read)])
  return(data.frame(meth = meth, depth = depth))
}, cl = 10) %>% bind_rows()

colnames(sample_output) <- paste0("sample_",colnames(sample_output))

write.table(bind_cols(marker_df, sample_output), file = "test_scripts/celfie_sample.txt", quote = F, row.names = F)

output_merged <- bind_cols(marker_df, sample_output, ref_merged)

#write.table(output_merged, file = "CelFiE_ref.txt", col.names = output_colnames, quote = F)
