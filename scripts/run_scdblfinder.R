suppressPackageStartupMessages({
  library(BiocParallel)
  library(Matrix)
  library(scDblFinder)
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: run_scdblfinder.R counts.mtx samples.tsv output.tsv seed expected_rate_or_auto")
}

counts <- readMM(args[[1]])
samples <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(counts) != nrow(samples)) {
  stop("Sample metadata rows do not match count-matrix columns")
}

set.seed(as.integer(args[[4]]))
sce <- SingleCellExperiment(assays = list(counts = counts))
colData(sce)$sample <- factor(samples$sample)

if (args[[5]] == "auto") {
  sce <- scDblFinder(sce, samples = colData(sce)$sample, BPPARAM = SerialParam())
} else {
  sce <- scDblFinder(
    sce,
    samples = colData(sce)$sample,
    dbr = as.numeric(args[[5]]),
    BPPARAM = SerialParam()
  )
}

result <- data.frame(
  barcode = samples$barcode,
  sample = samples$sample,
  score = colData(sce)$scDblFinder.score,
  class = as.character(colData(sce)$scDblFinder.class),
  stringsAsFactors = FALSE
)
write.table(result, args[[3]], sep = "\t", row.names = FALSE, quote = FALSE)

