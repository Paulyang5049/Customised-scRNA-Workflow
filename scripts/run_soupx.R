suppressPackageStartupMessages({
  library(Matrix)
  library(SoupX)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop(paste(
    "Usage: run_soupx.R raw.mtx filtered.mtx genes.txt raw_cells.txt",
    "cells.txt clusters.tsv corrected.mtx contamination.tsv seed"
  ))
}

raw_counts <- readMM(args[[1]])
filtered_counts <- readMM(args[[2]])
genes <- readLines(args[[3]])
raw_cells <- readLines(args[[4]])
cells <- readLines(args[[5]])
clusters <- read.delim(args[[6]], stringsAsFactors = FALSE, check.names = FALSE)

rownames(raw_counts) <- genes
colnames(raw_counts) <- raw_cells
rownames(filtered_counts) <- genes
colnames(filtered_counts) <- cells

set.seed(as.integer(args[[9]]))
soup_channel <- SoupChannel(raw_counts, filtered_counts, calcSoupProfile = TRUE)
cluster_vector <- setNames(clusters$cluster, clusters$barcode)
soup_channel <- setClusters(soup_channel, cluster_vector)
soup_channel <- autoEstCont(soup_channel, doPlot = FALSE)
corrected <- adjustCounts(soup_channel, roundToInt = TRUE)

writeMM(corrected, args[[7]])
contamination <- data.frame(
  barcode = rownames(soup_channel$metaData),
  contamination_fraction = soup_channel$metaData$rho,
  stringsAsFactors = FALSE
)
write.table(contamination, args[[8]], sep = "\t", row.names = FALSE, quote = FALSE)

