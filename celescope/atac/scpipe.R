library("scPipe")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--sample",help="Sample name", required=TRUE)
parser$add_argument("--r1",help="R1 fastq file", required=TRUE)
parser$add_argument("--r2",help="R2 fastq file", required=TRUE)
parser$add_argument("--r3",help="R3 fastq file", required=TRUE)
parser$add_argument("--thread",help="Sample name", required=TRUE)
parser$add_argument("--outdir", help='Output diretory.', required=TRUE)
parser$add_argument("--reference", help='Reference path', required=TRUE)
parser$add_argument("--organism", help='Organism', required=TRUE)
parser$add_argument("--rm_pcr_dup", help='Whether or not to remove PCR duplicate', required=TRUE)
args <- parser$parse_args()

sample <- args$sample
r1 <- args$r1 
r2 <- args$r2
r3 <- args$r3
thread <- as.double(args$thread)
output_folder <- args$outdir
reference <- args$reference
organism <- args$organism
rm_pcr_dup <- args$rm_pcr_dup

# convert format
sc_atac_trim_barcode (r1            = r1, 
                      r2            = r3, 
                      bc_file       = r2, 
                      rmN           = FALSE,
                      rmlow         = FALSE,
                      output_folder = output_folder)

# Aligning reads to a reference genome
# demux_r1 <- file.path(output_folder, paste0("demux_completematch_", sample, "_R1.fastq.gz"))
# demux_r2 <- file.path(output_folder, paste0("demux_completematch_", sample, "_R3.fastq.gz"))
demux_r1 <- file.path(output_folder, paste0("demux_completematch_", sample, "_1.fq"))
demux_r2 <- file.path(output_folder, paste0("demux_completematch_", sample, "_3.fq"))
aligned_bam <- sc_aligning(ref = reference,
                R1 = demux_r1,
                R2 = demux_r2,
                nthreads = thread,
                output_folder = output_folder)

# Demultiplexing the BAM file
sorted_tagged_bam <- sc_atac_bam_tagging(inbam = aligned_bam,
                       output_folder =  output_folder,
                       bam_tags = list(bc="CB", mb="OX"),
                       nthreads =  thread)

# Remove duplicates
removed <- sc_atac_remove_duplicates(sorted_tagged_bam,
                          output_folder = output_folder,)

if (grepl("True", rm_pcr_dup)){
  sorted_tagged_bam <- removed
}

# Gemerating a fragment file
sc_atac_create_fragments(inbam = sorted_tagged_bam,
                         output_folder = output_folder)

# Peak calling
sc_atac_peak_calling(inbam = sorted_tagged_bam,
                     ref = reference,
                     genome_size = NULL,
                     output_folder = output_folder)

# Assigning reads to features and feature counting
features          <- file.path(output_folder, "NA_peaks.narrowPeak")
sc_atac_feature_counting (fragment_file = file.path(output_folder, "fragments.bed"),
                          feature_input = features,
                          bam_tags      = list(bc="CB", mb="OX"),
                          feature_type  = "peak",
                          organism      = organism,
                          cell_calling  = "none",
                          min_uniq_frags = 0,
                          min_frac_peak = 0,
                          min_frac_promoter = 0,
                          yieldsize     = 1000000,
                          exclude_regions = TRUE,
                          output_folder = output_folder,
                          fix_chr       = "none",
                          create_report = FALSE
                          )

feature_matrix <- readRDS(file.path(output_folder, "unfiltered_feature_matrix.rds"))
dplyr::glimpse(feature_matrix)
sparseM <- readRDS(file.path(output_folder, "sparse_matrix.rds"))
dplyr::glimpse(sparseM)

# Generating the Single-cell Experiment (SCE) object
sce <- sc_atac_create_sce(input_folder = output_folder,
                   organism = organism,
                   feature_type = "peak",
                   pheno_data = NULL,
                   report = FALSE)