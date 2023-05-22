#!/usr/bin/env Rscript --vanilla


packagelist <- c("BiocManager",
                 "ggplot2",
                 "data.table",
                 "matrixStats",
                 "parallel",
                 "optparse",
                 "progress")

for(package in packagelist){
  if(!require(package, character.only = T, quietly = T)){ # check if package installed
    install.packages(package, quiet = T) # try installing from CRAN
    if(!require(package, character.only = T, quietly = T)){ # If package is still not installed
      BiocManager::install(package, ask = F) # try installing from Bioc
    }
    loaded_package <- require(package, character.only = T, quietly = T) # load package
  }
}

rm(list = ls())


# default values
def.bins.dir <- '/home/ssj36/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/train_test_split/split_train/'
def.markers.dir <- '/home/ssj36/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/train_test_split/split_train/markers/'
def.groups <- '/home/ssj36/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/train_test_split/split_train/groups/groups_test.csv'
def.blocks <- '/home/ssj36/Data/TestMixIn/from_megan_rad/Reference/ALL/CopyPAT/final_PAT_reference_SJ_Jan2023/train_test_split/split_train/train_blocks/train_blocks.bed.gz'
def.beta.thresh <- 0.3
def.bg.quant <- 0.10
def.target.quant <- 0.20


# parse command line arguments
parse.args <- function() {
  option_list = list(
    make_option("--bins_dir", type="character", default=def.bins.dir, 
                help="lbeta dir [default= %default]", metavar="character"),
    make_option(c("-o", "--out_dir"), type="character", default=def.markers.dir, 
                help="output dir [default= %default]", metavar="character"),
    make_option(c("-g", "--groups"), type="character", default=def.groups, 
                help="groups csv file [default= %default]", metavar="character"),
    make_option(c("-b", "--blocks"), type="character", default=def.blocks, 
                help="blocks file [default= %default]", metavar="character"),
    make_option(c("-t", "--target"), type="character", default='all', 
                help="group target, or all groups ('all') [default= %default]", metavar="character"),
    make_option('--top', type="integer", default=1000000, 
                help="output top TOP markers, sorted by margin [default= %default]"),
    make_option(c("-d", "--debug"), action="store_true", default=FALSE,
                help="Debug mode"),
    make_option(c("--no_header"), action="store_true", default=FALSE,
                help="output marker file header (e.g. #> Blood-B..."),
    make_option(c("--margin"), type="double", metavar="number", default=def.beta.thresh,
                help="beta difference threshold [default= %default]"),
    make_option(c("--bg_quant"), type="double", metavar="number", default=def.bg.quant,
                help="quantile relaxation of background group [default= %default]"),
    make_option(c("--target_quant"), type="double", metavar="number", default=def.target.quant,
                help="quantile relaxation of target group [default= %default]")
  ); 
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  opt
}
opt <- parse.args()
#blocks.path <- '/Users/nloyfer/work/atlas2/cell_paper/blocks.s150.p15.bed.gz'
blocks.path <- opt$blocks
bins.dir <- opt$bins_dir
markers.dir <- opt$out_dir
short.names.path <- opt$groups
req_target <- opt$target
# short.names.path <- 'all_bins/debug_all_short_names.csv'
min.cpgs <- 5
max.bp <- 2000
min.cov <- 10
beta.thresh <- opt$margin
target.quant <- opt$target_quant
bg.quant <- opt$bg_quant

output.header <- (!(opt$no_header))
DEBUG = opt$debug
nrows <- ifelse(DEBUG, 100000, -1)


# load blocks
load.blocks <- function(blocks.path, nrows = -1) {
  col.names <- c('chr', 'start', 'end', 'startCpG', 'endCpG') #,'anno', 'gene')
  df <- fread(blocks.path, col.names = col.names, nrows = nrows, fill = T)
}

load.sdf <- function() {
  # load short names table
  sdf <- data.table(read.table(short.names.path, sep = ',', header = T))
  if("include" %in% colnames(sdf)) {
    sdf <- sdf[toupper(include) == 'TRUE']
  }
  sdf
}

load.bins <- function(keep.inds) {
  sdf <- load.sdf()
  # init data.tables
  # res.beta - beta values table
  res.beta <- data.table(matrix(NA, nrow = sum(keep.inds), ncol = nrow(sdf)))
  colnames(res.beta) <- sdf[, name]
  # res.cov - coverage table
  res.cov <- copy(res.beta)
  
  # fill res.beta and res.cov column by column:
  cat(paste0('Loading ', nrow(sdf), ' bin files\n'))
  pb <- progress_bar$new(total = nrow(sdf))
  
  # row <- 1
  for (row in 1:nrow(sdf)) {
    pb$tick()
    name <- sdf[row, 'name'][[1]]
    rname <- paste0(gsub("\\+", ".", name), '.*lbeta')
    paths <- list.files(bins.dir, rname, full.names = T)
    if (length(paths) != 1) {
      print(paste0('Error matching name: ', name))
      quit(1)
    }
    path <- paths[1]
    #print(path)
    
    # load binary file to a data.table
    N <- file.info(path)$size / 2
    if (DEBUG) {
      N = nrows * 2
    }
    content <- data.table(matrix(readBin(con = path, what = "integer", 
                                         n = N, size = 2, signed=FALSE),
                                 nrow = N / 2, ncol = 2, byrow=TRUE))[keep.inds]
    
    colnames(content) <- c('meth', 'total')
    
    # update columns of res.beta and res.cov
    data.table::set(res.beta, j = which(name == colnames(res.beta)), value = content[, meth / total])
    data.table::set(res.cov, j = which(name == colnames(res.cov)), value = content[, total])
    rm(content)
  }
  list(res.beta, res.cov)
}


load.filtered.table <- function() {
  # load blocks
  df <- load.blocks(blocks.path, nrows)
  dim(df)
  # remove sex chromosomes
  keep.inds <- df$chr != 'chrM' # & df$chr != 'chrX' & df$chr!= 'chrY'
  
  # filter by block lengths
  keep.inds <- keep.inds & df$endCpG - df$startCpG >= min.cpgs & df$end - df$start <= max.bp
  
  
  # load binary files
  f <- load.bins(keep.inds)
  res.beta <- f[1][[1]]
  res.cov <- f[2][[1]]
  print(paste(nrow(df), "blocks before len filter"))
  df <- df[keep.inds]
  print(paste(nrow(df), "blocks left after len filter"))
  
  # Filter by coverage - average coverage of >= 10x
  avg.cov <- res.cov / (df$endCpG - df$startCpG)
  ss <- rowMeans(as.matrix(avg.cov)) >= min.cov
  res.beta <- res.beta[ss]
  res.cov <- res.cov[ss]
  df <- df[ss]
  print(paste(nrow(df), "blocks left after coverage filter"))
  
  # Filter by max-min difference
  #TODO - Revisit and test this filter
  ss <- (rowMaxs(as.matrix(res.beta, na.rm = T)) - 
           rowMins(as.matrix(res.beta, na.rm = T))) >= beta.thresh
  
  res.beta <- res.beta[ss]
  res.cov <- res.cov[ss]
  df <- df[ss]
  df$blockInd <- seq(nrow(df))
  print(paste(nrow(df), "blocks left after max-min filter"))
  
  list(df, res.beta, res.cov)
}


single.target <- function(target) {
  
  #vector of target samples
  target_samples <- sdf$group == target
  
  # filter by minimum coverage in target samples
  cov.table <- res.cov[, ..target_samples] / (dfo$endCpG - dfo$startCpG)
  ss <- rowMeans(as.matrix(cov.table)) >= min.cov
  df <- dfo[ss]
  res.beta.red <- res.beta[ss]
  
  
  # calc quantiles for target and background
  calc.quantiles <- function(table, rate) {
    q <- data.frame(rowQuantiles(as.matrix(
      table), na.rm = T, probs = c(rate, 1 - rate)))
    colnames(q) <- c('low', 'high')
    round(q, 2)
  }
  
  # quantile for target:
  #takes in filtered beta table for target and user specified target quantile
  qtarget <- calc.quantiles(res.beta.red[, ..target_samples], target.quant)
  
  # filter markers with no chance (e.g target$high > bg$max)
  ff <- (qtarget$high <= rowMaxs(as.matrix(res.beta.red[, !..target_samples]), na.rm = T) - beta.thresh) |
    (qtarget$low >= rowMins(as.matrix(res.beta.red[, !..target_samples]), na.rm = T) + beta.thresh)
  qtarget <- qtarget[ff, ]
  df <- df[ff]
  
  # quantile for background:
  qbg <- calc.quantiles(res.beta.red[ff, !..target_samples], bg.quant)
  
  # add quantiles info to df:
  df$target <- paste0(qtarget$low, '-', qtarget$high)
  df$bg <- paste0(qbg$low, '-', qbg$high)
  
  # find hyper/hypo markers
  hyper.cond <- qtarget$low - qbg$high >= beta.thresh
  hypo.cond <- qbg$low - qtarget$high >= beta.thresh
  
  # add direction (+/-) column
  df$direction <- NA
  df$direction[hyper.cond] <- '+'
  df$direction[hypo.cond] <- '-'
  
  # add margin column
  df$margin <- round(rowMaxs(cbind(qtarget$low - qbg$high, qbg$low - qtarget$high)), 2)
  #print(paste(nrow(df), "blocks left after quantile filter"))
  df <- na.omit(df)
  #print(paste(nrow(df), "blocks left after na.omit"))
  cat(paste0(nrow(df), ' markers\n'))
  if (nrow(df) == 0) {
    return()
  }
  
  # clip to top TOP markers
  #TODO - test this sort function, consider sorting by wilcoxon p-value
  if (nrow(df) > opt$top) {
    # sort by margin
    df <- df[order(-df$margin), ][1:opt$top, ]
    df <- df[order(df$startCpG), ]
  }
  
  # add blocks info
  # df$name <- paste0(df$chr, ":", df$start, '-', df$end, '.', target)
  df$name <- paste0(df$chr, ":", df$start, '-', df$end)
  df$label <- target
  df$CpGs <- paste0(df$endCpG - df$startCpG, 'CpGs')
  df$bp <- paste0(df$end - df$start, 'bp')
  
  names(df)[names(df) == 'regul'] <- 'anno'
  
  # add Wilcoxan p-value
  df$p.value <- NA
  if (sum(target_samples) == 1) {
    df$p.value <- NA
    } else {
      reduced.beta <- res.beta[df$blockInd]
      df$p.value <- mcmapply(function(i)
      {wilcoxon.test(reduced.beta[i, ..target_samples], reduced.beta[i, !..target_samples])$p.value}, 1:nrow(df))
      }
  df$q.value <- p.adjust(df$p.value, method = fdr, n = length(df$p.value))
  
  # add sign to margin and sort from + to - 
  #df$margin[df$direction == '-'] <- df$margin[df$direction == '-'] * (-1)
  #df <- df[order(df$margin), ]

  #Sort by adjusted p.value
  df <- df[order(df$q.value, decreasing = FALSE), ]
  
  # Dump
  dump.markers(df, sdf, target)
}

dump.markers <- function(df, sdf, target) {
  df <- df[, c('chr', 'start', 'end', 'startCpG', 'endCpG', 'label', 'name',
               'CpGs', 'bp', 'target', 'bg', 'margin', 'direction', #'anno', 
                'p.value', 'q.value')] #,'gene')]
  
  out.name <- file.path(markers.dir, paste0('markers.', target, '.bed'))
  cat(out.name, sep='\n')
  
  if (output.header) {
    # header - write target and bg sample name
    write.table(paste0('# >', sdf$name[sdf$group == target]), file=out.name, quote = F,row.names=F, col.names = FALSE)
    write.table(paste0('# <', sdf$name[sdf$group != target]), file=out.name, quote = F,row.names=F, col.names = FALSE, append = T)
  }
  
  fwrite(df, out.name, quote = F, na = 'NA', sep='\t', append = output.header)
}

# create markers dir, if not exists:
tmp <- ifelse(!dir.exists(markers.dir), dir.create(markers.dir), FALSE)

# load all data and metadata
r <- load.filtered.table()
res.beta <- r[2][[1]]
res.cov <- r[3][[1]]
dfo <- r[1][[1]]
sdf <- load.sdf()

groups <- sort(unique(sdf$group))
#req_target <- 'Endothel'

# single.target('Endothel')

for (cur_target in groups) {
  if (cur_target == 'bg') {
    next
  }
  if ((req_target != 'all') && (cur_target != req_target)) {
    next
  }
  cat(paste0('group ', cur_target, ':\t'), sep='\t')
  single.target(cur_target)  
  
}

