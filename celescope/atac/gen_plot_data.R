library(Seurat)
library(argparse)
library(Signac)
library(Gmisc)
library(ggplot2)

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
SeuratObj <- CreateSeuratObject(filtered_peak_count, project = sample, min.cells = 0, min.features = 0, assay = "ATAC")
SeuratObj <- subset(SeuratObj, nCount_ATAC != 0)

message("LSI analysis ...")
SeuratObj <- fastDoCall("RunTFIDF", c(object = SeuratObj))
SeuratObj <- FindTopFeatures(object = SeuratObj, min.cutoff = 'q0')
SeuratObj <- fastDoCall("RunSVD", c(object = SeuratObj))

message("UMAP analysis ...")
SeuratObj <- RunUMAP(object = SeuratObj, reduction = "lsi", dims = 1:30)
SeuratObj <- fastDoCall("FindNeighbors", c(object = SeuratObj, reduction = "lsi", dims = 1:30))
SeuratObj <- fastDoCall("FindClusters", c(object = SeuratObj, resolution = 0.6))
p1 <- DimPlot(object = SeuratObj, pt.size = 0.5, label = TRUE)
ggsave(file.path(outdir, paste0(sample, "_cluster.png")), p1, width=5, height=4)

df = merge(SeuratObj@reductions$umap@cell.embeddings, SeuratObj@meta.data, by = "row.names")
write.csv(df, meta_data, quote=FALSE, row.names=FALSE)                     
saveRDS(SeuratObj, rds)

