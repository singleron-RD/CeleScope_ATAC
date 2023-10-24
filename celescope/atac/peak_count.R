library(Seurat)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--peak_count_h5",help="peak count h5 file.", required=TRUE)
parser$add_argument("--peak_count_df",help="peak count txt file.", required=TRUE)
args <- parser$parse_args()
peak_count_h5 <- args$peak_count_h5 
peak_count_df <- args$peak_count_df

df <- Read10X_h5(peak_count_h5)
df <- as.data.frame(colSums(df))
df$barcode <- rownames(df)
df = df[,c(2,1)]
colnames(df) <- c("barcode", "peak_count")
write.csv(df, peak_count_df, quote=FALSE, row.names=FALSE)