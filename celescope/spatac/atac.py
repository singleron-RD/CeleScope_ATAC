import json
import subprocess
import pandas as pd
import numpy as np
import os
import collections
import tables
import scipy.sparse as sp_sparse
import scanpy as sc
import matplotlib.pyplot as plt
import snapatac2 as snap
from pybedtools import BedTool
from matplotlib import colors
from PIL import Image
from celescope.__init__ import ROOT_PATH
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.tools import get_plot_elements
from celescope.tools.plotly_plot import Insert_plot
from celescope.spatac.spatial import Spatial


__SUB_STEPS__ = ["mapping", "cells"]


FeatureBCMatrix = collections.namedtuple(
    "FeatureBCMatrix", ["ids", "names", "barcodes", "matrix"]
)

"""snapatac2 default variable name
BIN_SIZE = 500
N_FEATURES = 25000
N_COMPS = 30
N_NEIGHBORS = 30
MIN_DIST = 0.5
RESOLUTION = 0.5
MIN_COUNTS = 100
MIN_TSSE = 2
"""


def read_10X_h5(filename):
    """Read 10X HDF5 files, support both gene expression and peaks."""
    with tables.open_file(filename, "r") as f:
        try:
            group = f.get_node(f.root, "matrix")
        except tables.NoSuchNodeError:
            print("Matrix group does not exist in this file.")
            return None
        feature_group = getattr(group, "features")
        ids = getattr(feature_group, "id").read()
        names = getattr(feature_group, "name").read()
        barcodes = getattr(group, "barcodes").read()
        data = getattr(group, "data").read()
        indices = getattr(group, "indices").read()
        indptr = getattr(group, "indptr").read()
        shape = getattr(group, "shape").read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return FeatureBCMatrix(ids, names, barcodes, matrix)


def hires_nocrop_spatial(adata, **kwargs):
    """
    sc.pl.spatial wrapper: show full hires image (no crop)
    """
    spatial = adata.uns["spatial"]["sample1"]
    img = spatial["images"]["hires"]
    h, w = img.shape[:2]

    # scale factor (align coords)
    scale = spatial["scalefactors"].get("tissue_hires_scalef", 1.0)

    return sc.pl.spatial(
        adata,
        img_key="hires",
        show=False,
        crop_coord=(0, w / scale, 0, h / scale),
        **kwargs,
    )


def get_opts_atac(parser, sub_program):
    parser.add_argument(
        "--reference",
        help="Genome reference directory includes fasta, index, promoter bed file.",
        required=True,
    )
    parser.add_argument(
        "--genomesize",
        help="refer to www.genomesize.com. for example, 2.7e+9 for hs, 1.87e+9 for mm, 1.2e+8 for fruitfly",
        required=True,
    )
    parser.add_argument(
        "--peak_cutoff",
        help="Minimum number of peaks included in each cell",
        default="1",
    )
    parser.add_argument(
        "--count_cutoff",
        type=int,
        help="Cutoff for the number of count in each cell",
        default=1,
    )
    parser.add_argument(
        "--frip_cutoff",
        type=float,
        help="Cutoff for fraction of reads in promoter in each cell",
        default=0.01,
    )
    parser.add_argument(
        "--cell_cutoff",
        type=int,
        help="Minimum number of cells covered by each peak",
        default=1,
    )
    parser.add_argument(
        "--expected_target_cell_num",
        help="Expected cell number. If `--peak_cutoff` is provided, this argument is ignored.",
        type=int,
        default=8000,
    )
    parser.add_argument(
        "--coef",
        help="coef for auto peak cutoff, If `--peak_cutoff` is provided, this argument is ignored.",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--keep_mt", help="keep mitochondrial gene in fragments. ", action="store_true"
    )
    parser.add_argument("--out_sam", help="SAM file output. ", action="store_true")
    # snapatac2 variable name
    parser.add_argument(
        "--bin_size",
        help="The width, in base pairs, of each consecutive genomic bin.",
        type=int,
        default=500,
    )
    parser.add_argument(
        "--n_features",
        help="Maximum number of features to keep. The final number can be smaller if too few features pass filtering or have nonzero counts.",
        type=int,
        default=25000,
    )
    parser.add_argument(
        "--n_comps",
        help="Maximum number of spectral dimensions to compute.",
        type=int,
        default=30,
    )
    parser.add_argument(
        "--n_neighbors",
        help="Number of nearest neighbors to store for each observation.",
        type=int,
        default=30,
    )
    parser.add_argument(
        "--min_counts",
        help="Minimum data.obs.n_fragment value required for a cell to pass filtering.",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--min_tsse",
        help="Minimum data.obs.tsse value required for a cell to pass filtering.",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--min_dist",
        help="The min_dist parameter controls how tightly UMAP is allowed to pack points together.",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "--resolution",
        help="Clustering resolution. Larger values usually produce more clusters.",
        type=float,
        default=0.5,
    )
    if sub_program:
        s_common(parser)
        parser.add_argument("--spatial_dir", help="spatial directory", required=True)
        parser.add_argument(
            "--input_path", help="input_path from Barcode step.", required=True
        )
    return parser


class ATAC(Step):
    """
    - Run ATAC.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.input_path = os.path.abspath(args.input_path)
        self.reference = os.path.abspath(args.reference)
        self.genomesize = args.genomesize

        # 10X "/SGRNJ06/randd/USER/cjj/celedev/atac/MAESTRO/test/20240403_10X/737K-cratac-v1_rev.txt"
        # cut-off
        self.peak_cutoff = args.peak_cutoff
        self.count_cutoff = args.count_cutoff
        self.frip_cutoff = args.frip_cutoff
        self.cell_cutoff = args.cell_cutoff
        self.expected_target_cell_num = args.expected_target_cell_num
        self.coef = args.coef
        # self.in_tissue_barcodes = Spatial(
        #     self.args.spatial_dir
        # ).get_in_tissue_barcodes()
        # self.in_tissue_barcodes = [x.replace("_", "") for x in self.in_tissue_barcodes]

    @utils.add_log
    def mapping(self):
        """run chromap and process fragments"""
        cmd = (
            f"chromap --preset atac "
            f"-x {self.reference}/genome.index -r {self.reference}/genome.fa "
            f"-1 {self.input_path}/{self.sample}_S1_L001_R1_001.fastq -2 {self.input_path}/{self.sample}_S1_L001_R3_001.fastq "
            f"-b {self.input_path}/{self.sample}_S1_L001_R2_001.fastq "
            f"-o fragments_corrected_dedup_count.tsv -t {self.thread} "
        )
        subprocess.check_call(cmd, shell=True)
        df = pd.read_csv(
            "fragments_corrected_dedup_count.tsv",
            header=None,
            sep="\t",
            names=["chrom", "chromStart", "chromEnd", "barcode", "count"],
        )
        mapped_reads = sum(df["count"])
        self.add_metric(
            name="Confidently Unique mapped read pairs", value=mapped_reads, show=False
        )

        cmd = "grep -Ev 'chrM|MT' fragments_corrected_dedup_count.tsv > tmp.tsv"
        subprocess.check_call(cmd, shell=True)
        df = pd.read_csv(
            "tmp.tsv",
            header=None,
            sep="\t",
            names=["chrom", "chromStart", "chromEnd", "barcode", "count"],
        )
        mt_reads = mapped_reads - sum(df["count"])
        self.add_metric(name="Non-nuclear read pairs", value=mt_reads, show=False)
        if self.args.keep_mt:
            cmd = "rm tmp.tsv"
        else:
            cmd = "mv tmp.tsv fragments_corrected_dedup_count.tsv"
        subprocess.check_call(cmd, shell=True)

        cmd = (
            "bgzip -c fragments_corrected_dedup_count.tsv > fragments_corrected_dedup_count.tsv.gz;"
            "tabix -p bed fragments_corrected_dedup_count.tsv.gz"
        )
        subprocess.check_call(cmd, shell=True)

        if self.args.out_sam:
            cmd = (
                f"chromap --preset atac "
                f"-x {self.reference}/genome.index -r {self.reference}/genome.fa "
                f"-1 {self.input_path}/{self.sample}_S1_L001_R1_001.fastq -2 {self.input_path}/{self.sample}_S1_L001_R3_001.fastq "
                f"-b {self.input_path}/{self.sample}_S1_L001_R2_001.fastq "
                f"--SAM -o alignment.sam -t {self.thread} "
            )
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def qcstat(self):
        """qc stat"""
        # qc mapping
        cmd1 = (
            "sort -k4,4 -V fragments_corrected_dedup_count.tsv > fragments_corrected_count_sortedbybarcode.tsv;"
            "bedtools groupby -i fragments_corrected_count_sortedbybarcode.tsv -g 4 -c 5,5 -o count,sum > singlecell_mapped.txt"
        )

        # qc promoter
        cmd2 = (
            f"bedtools intersect -wa -a fragments_corrected_dedup_count.tsv -b {self.reference}/promoter.bed -u > fragments_promoter.tsv;"
            "sort -k4,4 -V fragments_promoter.tsv > fragments_promoter_sortbybarcode.tsv;"
            "bedtools groupby -i  fragments_promoter_sortbybarcode.tsv -g 4 -c 5 -o sum > singlecell_promoter.txt"
        )

        # qc singlecell
        cmd3 = (
            "sort -k1,1 singlecell_mapped.txt > singlecell_mapped_sortbybarcode.txt;"
            "sort -k1,1 singlecell_promoter.txt > singlecell_promoter_sortbybarcode.txt;"
            "join --nocheck-order -t $'\t' -a1 -e'0' -o'1.1 1.2 1.3 2.2' -1 1 -2 1 \
                singlecell_mapped_sortbybarcode.txt singlecell_promoter_sortbybarcode.txt > singlecell.txt"
        )

        for cmd in [cmd1, cmd2, cmd3]:
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def call_peak(self):
        """peak calling"""
        cmd1 = (
            f"macs2 callpeak -f BEDPE -g {self.genomesize} --outdir peak -n {self.sample} "
            "-B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t fragments_corrected_dedup_count.tsv"
        )

        # used to call shortpeak
        cmd2 = (
            "awk -F'\\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if (abs($3-$2)<=150) print}}' \
                fragments_corrected_dedup_count.tsv > fragments_corrected_150bp.tsv"
        )

        # call shortpeak
        cmd3 = (
            f"macs2 callpeak -f BEDPE -g {self.genomesize} --outdir peak -n {self.sample}_150bp "
            "-B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t fragments_corrected_150bp.tsv"
        )

        # merge peak
        cmd4 = (
            f"cat peak/{self.sample}_peaks.narrowPeak peak/{self.sample}_150bp_peaks.narrowPeak \
                | sort -k1,1 -k2,2n | cut -f 1-4 > tmp_peaks.bed;"
            f"mergeBed -i tmp_peaks.bed > peak/{self.sample}_final_peaks.bed;"
            "rm tmp_peaks.bed"
        )

        for cmd in [cmd1, cmd2, cmd3, cmd4]:
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_valid_cells(self):
        """filter by count-cutoff and frip-cutfoff"""
        df = pd.read_csv(
            "singlecell.txt",
            header=None,
            sep="\t",
            names=[
                "barcode",
                "fragments",
                "readSupport",
                "fragments_overlapping_promoter",
            ],
        )
        df = df[df["readSupport"] >= self.count_cutoff]
        df["fraction_in_promoter"] = (
            df["fragments_overlapping_promoter"] / df["readSupport"]
        )
        df = df[df["fraction_in_promoter"] >= self.frip_cutoff]
        # df = df[df["barcode"].isin(self.in_tissue_barcodes)]
        df["barcode"].to_csv("validcells.txt", header=None, index=None)

    @utils.add_log
    def gen_peak_matrix(self):
        """generate peak count matrix h5"""

        cmd = (
            f"python {ROOT_PATH}/atac/scatac_peakcount.py --barcode validcells.txt "
            f"--count_cutoff {self.count_cutoff} --cores {self.thread} --outprefix {self.sample} "
            f"--fragment fragments_corrected_dedup_count.tsv --peak peak/{self.sample}_final_peaks.bed "
        )

        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def filter_peak_matrix(self):
        """filter peak count matrix h5"""
        raw_matrix = f"peak/{self.sample}_peak_count.h5"
        peak_cutoff = ATAC.count_peak_cutoff(
            raw_matrix, self.peak_cutoff, self.expected_target_cell_num, self.coef
        )
        print(f"Peak Cutoff: {peak_cutoff}")

        cmd = (
            f"python {ROOT_PATH}/atac/scatac_filter.py --peakcount {raw_matrix} "
            f"--peak_cutoff {peak_cutoff} --cell_cutoff {self.cell_cutoff} --outprefix {self.sample}"
        )

        subprocess.check_call(cmd, shell=True)

    @staticmethod
    def count_peak_cutoff(raw_matrix, peak_cutoff, expected_target_cell_num, coef):
        if peak_cutoff == "auto":
            scatac_count = read_10X_h5(raw_matrix)
            peakmatrix = scatac_count.matrix
            features = scatac_count.names.tolist()
            barcodes = scatac_count.barcodes.tolist()
            if isinstance(features[0], bytes):
                features = [i.decode() for i in features]
            if isinstance(barcodes[0], bytes):
                barcodes = [i.decode() for i in barcodes]

            peaks_per_cell = np.asarray((peakmatrix > 0).sum(axis=0))
            sorted_counts = sorted(peaks_per_cell[0], reverse=True)
            return int(
                np.percentile(sorted_counts[:expected_target_cell_num], 99) * coef
            )

        else:
            return int(peak_cutoff)

    def run(self):
        cwd = os.getcwd()
        os.chdir(self.outdir)
        self.mapping()
        self.qcstat()
        self.call_peak()
        self.get_valid_cells()
        self.gen_peak_matrix()
        self.filter_peak_matrix()
        os.chdir(cwd)


def atac(args):
    with ATAC(args) as runner:
        runner.run()

    with Maestro_metrics(args) as runner:
        runner.run()

    with Mapping(args) as runner:
        runner.run()

    with Cells(args) as runner:
        runner.run()


class Maestro_metrics(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.df_mapping = pd.read_csv(
            f"{self.outdir}/fragments_corrected_count_sortedbybarcode.tsv",
            header=None,
            sep="\t",
            names=["chrom", "chromStart", "chromEnd", "barcode", "count"],
        )
        self.df_barcode = pd.read_csv(
            f"{self.outdir}/singlecell.txt",
            header=None,
            sep="\t",
            names=["barcode", "fragments", "readSupport", "overlap_promoter"],
        )
        self.df_barcode["frac_promoter"] = round(
            self.df_barcode["overlap_promoter"] / self.df_barcode["readSupport"], 4
        )

        self.df_fragments = (
            f"{self.outdir}/fragments_corrected_count_sortedbybarcode.tsv"
        )
        self.df_peaks = f"{self.outdir}/peak/{self.sample}_final_peaks.bed"

        # out
        self.spatial = f"{self.outdir}/spatial"
        self.df_cell_metrics = f"{self.outdir}/cell_qc_metrics.tsv"
        self.df_umap_file = f"{self.outdir}/umap_coord.tsv"
        self.raw_counts_png = f"{self.outdir}/raw_counts.png"
        self.filtered_counts_png = f"{self.outdir}/counts.png"
        self.cluster_png = f"{self.outdir}/cluster.png"

    @utils.add_log
    def read_data(self):
        genome = snap.genome.Genome(
            fasta=f"{self.args.reference}/genome.fa",
            annotation=f"{self.args.reference}/gene.gtf",
        )
        fragment = f"{self.outdir}/fragments_corrected_dedup_count.tsv"
        self.adata = snap.pp.import_data(
            fragment,
            chrom_sizes=genome,
            min_num_fragments=0,
            sorted_by_barcode=False,
            chrM=["chrM", "M", "MT", "Mt"],
        )

        # QC
        snap.metrics.tsse(self.adata, genome)
        snap.pp.filter_cells(
            self.adata, min_counts=self.args.min_counts, min_tsse=self.args.min_tsse
        )

    @utils.add_log
    def run_snapatac2(self):
        snap.pp.add_tile_matrix(self.adata, bin_size=self.args.bin_size)
        snap.pp.select_features(self.adata, n_features=self.args.n_features)
        snap.tl.spectral(self.adata, n_comps=self.args.n_comps)
        snap.pp.knn(self.adata, n_neighbors=self.args.n_neighbors)
        snap.tl.leiden(self.adata, resolution=self.args.resolution)
        snap.tl.umap(
            self.adata, n_neighbors=self.args.n_neighbors, min_dist=self.args.min_dist
        )

    @utils.add_log
    def write_umap(self):
        df_umap = pd.DataFrame(
            self.adata.obsm["X_umap"],
            index=self.adata.obs_names,
            columns=[f"UMAP_{i+1}" for i in range(self.adata.obsm["X_umap"].shape[1])],
        )
        df_umap["cluster"] = self.adata.obs.leiden
        df_umap["cluster"] = df_umap["cluster"].map(lambda x: int(x) + 1)
        df_umap.to_csv(self.df_umap_file, sep="\t")

    @utils.add_log
    def write_spatial(self):
        coords = pd.read_csv(f"{self.spatial}/tissue_positions_list.csv", header=None)
        coords.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row",
            "pxl_col",
        ]
        coords["barcode"] = coords["barcode"].str.replace("_", "")
        coords = coords.set_index("barcode")
        common = self.adata.obs_names.intersection(coords.index)
        self.adata = self.adata[common]
        coords = coords.loc[common]

        # 写入空间坐标
        self.adata.obs["in_tissue"] = coords["in_tissue"].astype(int)
        self.adata.obsm["spatial"] = coords[["pxl_col", "pxl_row"]].values

    @utils.add_log
    def add_fragment_count(self):
        # self.adata.obs["fragment_counts"] = list(self.df_barcode["fragments"])
        self.adata.obs["fragment_counts"] = self.adata.obs.index.map(
            self.df_barcode.set_index("barcode")["fragments"]
        )
        self.adata.obs["log_fragments"] = np.log1p(self.adata.obs["fragment_counts"])

    @utils.add_log
    def add_visium(self):
        hires = Image.open(f"{self.spatial}/tissue_hires_image.png").convert("RGB")
        lowres = Image.open(f"{self.spatial}/tissue_lowres_image.png").convert("RGB")

        hires = np.array(hires)
        lowres = np.array(lowres)

        with open(f"{self.spatial}/scalefactors_json.json") as f:
            scalefactors = json.load(f)

        library_id = "sample1"

        self.adata.uns["spatial"] = {
            library_id: {
                "images": {"hires": hires, "lowres": lowres},
                "scalefactors": scalefactors,
                "metadata": {},
            }
        }

    @utils.add_log
    def add_count_plot(self, plot_path):
        plt.figure(figsize=(8, 8))
        hires_nocrop_spatial(
            self.adata,
            color=["fragment_counts"],
            color_map="Reds",
            size=1.5,
            alpha=0.5,
            norm=colors.LogNorm(vmin=1),
        )
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()

    @utils.add_log
    def add_cluster_plot(self, plot_path):
        plt.figure(figsize=(8, 8))
        hires_nocrop_spatial(self.adata, color=["leiden"], size=1.5)
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()

    @utils.add_log
    def get_filtered_data(self):
        adata_filtered = self.adata[self.adata.obs["in_tissue"] == 1, :].copy()
        self.adata = adata_filtered

    @utils.add_log
    def run(self):
        self.read_data()
        spatial = Spatial(self.args.spatial_dir)
        spatial.output_spatial(self.spatial)
        self.write_spatial()
        self.add_fragment_count()
        self.add_visium()
        self.add_count_plot(self.raw_counts_png)
        self.get_filtered_data()
        self.run_snapatac2()
        self.write_umap()
        self.add_count_plot(self.filtered_counts_png)
        self.add_cluster_plot(self.cluster_png)


class Mapping(Maestro_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.cell_barcode = utils.read_one_col(self.df_umap_file)[0]
        self.cell_barcode = [item.split("\t")[0] for item in self.cell_barcode]
        del self.cell_barcode[0]
        self.df_mapping_cell = self.df_mapping[
            self.df_mapping["barcode"].isin(self.cell_barcode)
        ]

    def run(self):
        valid_reads = self.get_slot_key(
            slot="metrics",
            step_name="barcode",
            key="Valid Reads",
        )
        mapped_reads = self.get_slot_key(
            slot="metrics",
            step_name="aTAC",
            key="Confidently Unique Mapped Read Pairs",
        )
        mt_reads = self.get_slot_key(
            slot="metrics",
            step_name="aTAC",
            key="Non-nuclear Read Pairs",
        )

        self.add_metric(
            name="Confidently Unique mapped read pairs",
            value=mapped_reads,
            total=valid_reads,
            help_info="Fraction of Unique mapped read pairs",
        )

        self.add_metric(
            name="Non-nuclear read pairs",
            value=mt_reads,
            total=valid_reads,
            help_info="Fraction of Unique mapped read pairs that map to non-nuclear genome contigs, including mitochondria",
        )

        self.df_mapping["size"] = (
            self.df_mapping["chromEnd"] - self.df_mapping["chromStart"]
        )

        self.add_metric(
            name="Fragments in nucleosome-free regions",
            value=self.df_mapping[self.df_mapping["size"] <= 124].shape[0],
            total=self.df_mapping.shape[0],
            help_info="Fraction of high-quality fragments smaller than 124 basepairs.",
        )

        self.add_metric(
            name="Fragments flanking a single nucleosome",
            value=self.df_mapping[
                (self.df_mapping["size"] > 124) & (self.df_mapping["size"] <= 296)
            ].shape[0],
            total=self.df_mapping.shape[0],
            help_info="Fraction of high-quality fragments between 124 and 296 basepairs.",
        )

        self.add_metric(
            name="Percent duplicates",
            value=self.df_mapping_cell[self.df_mapping_cell["count"] > 1].shape[0],
            total=self.df_mapping_cell.shape[0],
            help_info="Fraction of high-quality read pairs that are deemed to be PCR duplicates. A high-quality read-pair is one with mapping quality > 30, that is not chimeric and maps to nuclear contigs. This metric is a measure of sequencing saturation and is a function of library complexity and sequencing depth. More specifically, this is the fraction of high-quality fragments with a valid barcode that align to the same genomic position as another read pair in the library.",
        )

        self.df_mapping = self.df_mapping[self.df_mapping["size"] < 800]
        Insertplot = Insert_plot(df=self.df_mapping).get_plotly_div()
        self.add_data(Insert_plot=Insertplot)


class Cells(Maestro_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.cell_barcode = utils.read_one_col(self.df_umap_file)[0]
        self.cell_barcode = [item.split("\t")[0] for item in self.cell_barcode]
        del self.cell_barcode[0]

    @utils.add_log
    def count_overlap_peak(self):
        """count fragments overlapping peaks"""
        fragments = BedTool(self.df_fragments)
        peaks = BedTool(self.df_peaks)
        frip = fragments.intersect(peaks, u=True, wa=True)
        df_in_peaks = frip.to_dataframe()

        final_df_count = df_in_peaks.groupby("name", as_index=False).agg(
            {"score": "count"}
        )
        final_df_count = final_df_count.rename(
            columns={"name": "barcode", "score": "count"}
        )
        self.df_barcode = pd.merge(
            self.df_barcode, final_df_count, on="barcode", how="outer"
        ).fillna(0)
        self.df_barcode = self.df_barcode.rename(columns={"count": "overlap_peaks"})
        self.df_barcode["frac_peak"] = round(
            self.df_barcode["overlap_peaks"] / self.df_barcode["fragments"], 4
        )
        self.df_barcode["cell_called"] = self.df_barcode["barcode"].apply(
            lambda x: True if x in self.cell_barcode else False
        )
        self.df_barcode.to_csv(self.df_cell_metrics, sep="\t", index=False)

    def run(self):
        self.count_overlap_peak()

        self.df_cell_barcode = self.df_barcode[
            self.df_barcode["barcode"].isin(self.cell_barcode)
        ]
        cell_num = len(self.cell_barcode)

        self.add_metric(
            name="In Tissue Spots",
            value=len(self.cell_barcode),
            help_info="Number of spots in tissue determined based on the image.",
        )

        raw_reads = self.get_slot_key(
            slot="metrics",
            step_name="barcode",
            key="Raw Reads",
        )
        self.add_metric(
            name="Mean raw read pairs per Spot",
            value=int(raw_reads / cell_num),
            help_info="Total number of read pairs divided by the number of in tissue spot.",
        )

        total_fragments = sum(self.df_barcode.readSupport)
        cell_fragments = sum(self.df_cell_barcode.readSupport)
        self.add_metric(
            name="Fraction of high-quality fragments in Spots",
            value=f"{round(cell_fragments / total_fragments * 100, 2)}%",
            help_info="Fraction of high-quality fragments with a in tissue spot that are associated with cell-containing partitions.",
        )

        frac_peak = self.df_cell_barcode["frac_peak"].mean()
        self.add_metric(
            name="Fraction of Fragments Overlap with Peaks in Spots",
            value=f"{round(frac_peak * 100, 2)}%",
            help_info="The proportion of fragments in a in tissue spot to overlap with a peak.",
        )

        frac_promoter = self.df_cell_barcode["frac_promoter"].mean()
        self.add_metric(
            name="Fraction of fragments overlap with promoter in Spots",
            value=f"{round(frac_promoter * 100, 2)}%",
            help_info="The proportion of fragments in a in tissue spot to overlap with a promoter sequence.",
        )

        self.add_metric(
            name="Median high-quality fragments per Spot",
            value=int(np.median(self.df_cell_barcode.fragments)),
            help_info="The median number of high-quality fragments per spot",
        )

        self.df_barcode.sort_values(by="overlap_peaks", ascending=False, inplace=True)
        self.df_barcode.iloc[0, -1] = True
        self.add_data(
            chart=get_plot_elements.plot_barcode_rank(
                self.df_barcode, log_uniform=False
            )
        )
