library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(argparse)
library(Seurat)

parser <- ArgumentParser()
parser$add_argument("--sce_rds",help="scPipe atac SCEobject rds file", required=TRUE)
parser$add_argument("--meta_data", help='meta-data file.', required=TRUE)
args <- parser$parse_args()

sce_rds <- args$sce_rds
meta_data <- args$meta_data

rds = readRDS(sce_rds)
df = merge(rds$ATAC@reductions$umap@cell.embeddings, rds$ATAC@meta.data, by = "row.names")
write.csv(df, meta_data, quote=FALSE, row.names=FALSE)
