library(MAESTRO)
library(Seurat)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--filtered_peak_count", help="filter peak count", required=TRUE)
parser$add_argument("--rds", help="rds file", required=TRUE)
parser$add_argument("--meta_data", help="meta data file", required=TRUE)
parser$add_argument("--outdir", help="output directory", required=TRUE)
parser$add_argument("--sample", help="sample name", required=TRUE)
args <- parser$parse_args()

filtered_peak_count <- args$filtered_peak_count
rds <- args$rds
meta_data <- args$meta_data
outdir <- args$outdir
sample <- args$sample

filtered_peak_count = Read10X_h5(filtered_peak_count)
rds.res <- ATACRunSeurat(inputMat = filtered_peak_count,
                                 project = sample,
                                 min.c = 0,
                                 min.p = 0,
                                 method = "LSI",
                                 dims.use = 1:30,
                                 cluster.res = 0.6,
                                 only.pos = TRUE,
                                 peaks.cutoff = 0,
                                 peaks.pct = 0,
                                 peaks.logfc = 0,
                                 outdir = outdir)
                                 
saveRDS(rds.res, rds)
df = merge(rds.res$ATAC@reductions$umap@cell.embeddings, rds.res$ATAC@meta.data, by = "row.names")
write.csv(df, meta_data, quote=FALSE, row.names=FALSE)
