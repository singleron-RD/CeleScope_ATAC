library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(argparse)
library(irlba)
library(RANN)
library(umap)
library(igraph)
library(tibble)
library(Rtsne)

parser <- ArgumentParser()
parser$add_argument("--feature_matrix",help="unfiltered feature matrix rds file.", required=TRUE)
parser$add_argument("--fragments",help="fragments bed file.", required=TRUE)
parser$add_argument("--insert_size_data",help="insert size data csv file", required=TRUE)
parser$add_argument("--binary_matrix",help="binary matrix rds file.", required=TRUE)
parser$add_argument("--sce_rds",help="scPipe atac SCEobject rds file", required=TRUE)
parser$add_argument("--cell_tsne_info", help='cell qc metrics with tsne coord csv file.', required=TRUE)
args <- parser$parse_args()

feature_matrix <- args$feature_matrix
fragments <- args$fragments 
insert_size_data <- args$insert_size_data
binary_matrix <- args$binary_matrix
sce_rds <- args$sce_rds
cell_tsne_info <- args$cell_tsne_info

# Insert size distribution
unfiltered_mtx <- readRDS(feature_matrix)
unfiltered_mtx_bcs <- colnames(unfiltered_mtx)
frags <- fread(fragments)[V4 %in% unfiltered_mtx_bcs, ]

frags[, 'isize' := V3 - V2]
if (nrow(frags) >= 100000) {
    frags = frags[sort(sample(1:nrow(frags), 100000)), ]
}
frags <- frags[isize < 800]
write.csv(frags, insert_size_data, quote=FALSE, row.names=FALSE)


TF.IDF.custom <- function(binary.mat, verbose = TRUE) {
    object <- binary.mat
    npeaks       <- Matrix::colSums(x = object)
    tf           <- Matrix::tcrossprod(x = as.matrix(object), y = Matrix::Diagonal(x = 1 / npeaks))
    rsums        <- Matrix::rowSums(x = object)
    idf          <- ncol(x = object) / rsums
    norm.data    <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
    scale.factor <- 1e4
    slot(object = norm.data, name = "x") <- log1p(x = slot(object = norm.data, name = "x") * scale.factor)
    norm.data[which(x = is.na(x = norm.data))] <- 0
    return(norm.data)
  }

binary.mat <- readRDS(binary_matrix)
sce <-readRDS(sce_rds)
counts <- assay(sce)
bin_mat <- as.matrix((counts>0)+0)
binary.mat <- TF.IDF.custom(bin_mat)

set.seed(123)
n_bcs <- max(min(50, ncol(binary.mat), nrow(binary.mat))-1,0)
mat.lsi          <- irlba(binary.mat, n_bcs)
d_diagtsne       <- matrix(0, n_bcs, n_bcs)
diag(d_diagtsne) <- mat.lsi$d
mat_pcs          <- t(d_diagtsne %*% t(mat.lsi$v))
# binary.mat -> bin_mat
rownames(mat_pcs)<- colnames(bin_mat)

# clustering in the PCA space using KNN --------------

library(RANN)
knn.info<- RANN::nn2(mat_pcs, k = 30)

## convert to adjacency matrix
knn           <- knn.info$nn.idx
adj           <- matrix(0, nrow(mat_pcs), nrow(mat_pcs))
rownames(adj) <- colnames(adj) <- rownames(mat_pcs)
for(i in seq_len(nrow(mat_pcs))) {
  adj[i,rownames(mat_pcs)[knn[i,]]] <- 1
}

## convert to graph
library(igraph)
g <- igraph::graph.adjacency(adj, mode="undirected")
g <- simplify(g) ## remove self loops

# identify communities, many algorithums. Use the Louvain clustering ------------
km         <- igraph::cluster_louvain(g)
com        <- km$membership
names(com) <- km$names

# # running UMAP ------------------------------
# norm.data.umap    <- umap::umap(mat_pcs)
# df_umap           <- as.data.frame(norm.data.umap$layout)
# colnames(df_umap) <- c("UMAP1", "UMAP2")
# df_umap$barcode   <- rownames(mat_pcs)

# running tSNE ------------------------------
norm.data.tsne    <- Rtsne(mat_pcs)
df_tsne           <- as.data.frame(norm.data.tsne$Y)
colnames(df_tsne) <- c("tSNE_1", "tSNE_2")
df_tsne$barcode   <- rownames(mat_pcs)

df_tsne <- dplyr::left_join(df_tsne, enframe(com), by = c("barcode" = "name")) %>%
  dplyr::rename(cluster = value) %>%
  dplyr::mutate(cluster = as.factor(cluster))

sce_coldata <- colData(sce)
  df_tsne <- base::merge(df_tsne, sce_coldata, by.x = "barcode", by.y = "row.names", all.x = TRUE) 

write.csv(df_tsne, cell_tsne_info, quote=FALSE, row.names=FALSE)
