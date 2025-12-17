import subprocess
import pandas as pd
import numpy as np
import os
import itertools
import collections
import tables
import scipy.sparse as sp_sparse
import scanpy as sc
from sklearn.decomposition import TruncatedSVD
from multiprocessing import Pool
from celescope.__init__ import ROOT_PATH
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.tools import get_plot_elements
from celescope.tools.plotly_plot import Insert_plot


__SUB_STEPS__ = ["mapping", "cells"]


FeatureBCMatrix = collections.namedtuple(
    "FeatureBCMatrix", ["ids", "names", "barcodes", "matrix"]
)


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
        default="auto",
    )
    parser.add_argument(
        "--count_cutoff",
        type=int,
        help="Cutoff for the number of count in each cell",
        default=1000,
    )
    parser.add_argument(
        "--frip_cutoff",
        type=float,
        help="Cutoff for fraction of reads in promoter in each cell",
        default=0.2,
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
    if sub_program:
        s_common(parser)
        parser.add_argument(
            "--match_dir", help="Matched scRNA-seq directory", required=True
        )
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
        self.match_dir = args.match_dir

        # 10X "/SGRNJ06/randd/USER/cjj/celedev/atac/MAESTRO/test/20240403_10X/737K-cratac-v1_rev.txt"
        # cut-off
        self.peak_cutoff = args.peak_cutoff
        self.count_cutoff = args.count_cutoff
        self.frip_cutoff = args.frip_cutoff
        self.cell_cutoff = args.cell_cutoff
        self.expected_target_cell_num = args.expected_target_cell_num
        self.coef = args.coef

        if self.match_dir != "None":
            self.match_cell_barcodes, _ = utils.get_barcode_from_match_dir(
                self.match_dir
            )
            self.match_cell_barcodes = [
                item.replace("_", "") for item in self.match_cell_barcodes
            ]
            self.peak_cutoff = 0
            self.count_cutoff = 0
            self.frip_cutoff = 0
            self.cell_cutoff = 0

    @utils.add_log
    def mapping(self):
        """run chromap and process fragments"""
        cmd1 = (
            f"chromap --preset atac "
            f"-x {self.reference}/genome.index -r {self.reference}/genome.fa "
            f"-1 {self.input_path}/{self.sample}_S1_L001_R1_001.fastq -2 {self.input_path}/{self.sample}_S1_L001_R3_001.fastq "
            f"-b {self.input_path}/{self.sample}_S1_L001_R2_001.fastq "
            f"-o fragments_corrected_dedup_count.tsv -t {self.thread} "
        )

        cmd2 = "grep -Ev 'chrM|MT' fragments_corrected_dedup_count.tsv > tmp && mv tmp fragments_corrected_dedup_count.tsv"

        cmd3 = (
            "bgzip -c fragments_corrected_dedup_count.tsv > fragments_corrected_dedup_count.tsv.gz;"
            "tabix -p bed fragments_corrected_dedup_count.tsv.gz"
        )

        cmds = [cmd1, cmd2, cmd3]
        if self.args.keep_mt:
            del cmds[1]

        for cmd in cmds:
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def qcstat(self):
        """qc stat"""
        # qc mapping
        cmd1 = (
            "sort -k4,4 -V fragments_corrected_dedup_count.tsv > fragments_corrected_count_sortedbybarcode.tsv;"
            "bedtools groupby -i fragments_corrected_count_sortedbybarcode.tsv -g 4 -c 5 -o sum > singlecell_mapped.txt"
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
            "join --nocheck-order -t $'\t' -a1 -e'0' -o'1.1 1.2 2.2' -1 1 -2 1 \
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
            names=["barcode", "fragments", "fragments_overlapping_promoter"],
        )
        df = df[df["fragments"] >= self.count_cutoff]
        df["fraction_in_promoter"] = (
            df["fragments_overlapping_promoter"] / df["fragments"]
        )
        df = df[df["fraction_in_promoter"] >= self.frip_cutoff]
        if self.match_dir != "None":
            df = df[df["barcode"].isin(self.match_cell_barcodes)]
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

        self.filtered_peak_count = (
            f"{self.outdir}/peak/{self.sample}_filtered_peak_count.h5"
        )
        self.adata = sc.read_10x_h5(self.filtered_peak_count, gex_only=False)
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
            names=["barcode", "fragments", "overlap_promoter"],
        )
        self.df_barcode["frac_promoter"] = round(
            self.df_barcode["overlap_promoter"] / self.df_barcode["fragments"], 4
        )

        self.df_fragments = pd.read_csv(
            f"{self.outdir}/fragments_corrected_count_sortedbybarcode.tsv",
            header=None,
            sep="\t",
            names=["chr", "start", "end", "barcode", "count"],
        )
        self.df_peaks = pd.read_csv(
            f"{self.outdir}/peak/{self.sample}_final_peaks.bed",
            header=None,
            sep="\t",
            names=["chr", "start", "end"],
        )

        # out
        self.df_cell_metrics = f"{self.outdir}/cell_qc_metrics.tsv"
        self.df_tsne_file = f"{self.outdir}/tsne_coord.tsv"

    @utils.add_log
    def run_tfidf(self):
        """Running TF-IDF normalization"""
        layer_name = "tfidf"
        X = self.adata.X
        if not sp_sparse.issparse(X):
            X = sp_sparse.csr_matrix(X)

        # TF-IDF
        cell_sum = np.array(X.sum(axis=1)).flatten()
        tf = X.multiply(1 / cell_sum[:, None])
        N = X.shape[0]
        peak_n = np.array((X > 0).sum(axis=0)).flatten()
        idf = np.log1p(N / (1 + peak_n))
        tfidf = tf.multiply(idf)
        self.adata.layers[layer_name] = tfidf.tocsr()
        # log1p
        sc.pp.log1p(self.adata)

    @utils.add_log
    def run_hvg(self):
        """sc.highly_variable_genes()"""
        sc.pp.highly_variable_genes(self.adata)

    @utils.add_log
    def run_svd(self):
        """Running SVD"""
        n_components = 30
        layer = "tfidf"
        key_added = "X_lsi"
        X = self.adata.layers[layer]
        if not sp_sparse.issparse(X):
            X = sp_sparse.csr_matrix(X)

        svd = TruncatedSVD(
            n_components=n_components, random_state=0, algorithm="randomized", n_iter=3
        )
        lsi = svd.fit_transform(X)
        self.adata.obsm[key_added] = lsi
        self.adata.uns["lsi"] = {"svd": svd}

    @utils.add_log
    def run_neighbors(self):
        """Computing neighbors using LSI"""
        sc.pp.neighbors(self.adata, use_rep="X_lsi", n_neighbors=30, n_pcs=30)

    @utils.add_log
    def run_tsne(self):
        """Running TSNE"""
        sc.tl.tsne(self.adata, use_rep="X_lsi", n_pcs=30, perplexity=30, random_state=0)

    @utils.add_log
    def run_leiden(self):
        """FindClusters (Leiden)"""
        sc.tl.leiden(self.adata, resolution=0.8)

    @utils.add_log
    def write_tsne(self):
        df_tsne = self.adata.obsm.to_df()[["X_tsne1", "X_tsne2"]]
        df_tsne["cluster"] = self.adata.obs.leiden
        tsne_name_dict = {"X_tsne1": "tSNE_1", "X_tsne2": "tSNE_2"}
        df_tsne = df_tsne.rename(tsne_name_dict, axis="columns")
        df_tsne["cluster"] = df_tsne["cluster"].map(lambda x: int(x) + 1)
        df_tsne.to_csv(self.df_tsne_file, sep="\t")

    @utils.add_log
    def run(self):
        self.run_tfidf()
        self.run_hvg()
        self.run_svd()
        self.run_neighbors()
        self.run_tsne()
        self.run_leiden()
        self.write_tsne()


class Mapping(Maestro_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.cell_barcode = utils.read_one_col(self.df_tsne_file)[0]
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

        self.df_mapping["size"] = (
            self.df_mapping["chromEnd"] - self.df_mapping["chromStart"]
        )

        self.add_metric(
            name="Confidently Unique mapped read pairs",
            value=sum(self.df_mapping["count"]),
            total=valid_reads,
            help_info="Fraction of Unique mapped read pairs",
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

        self.cell_barcode = utils.read_one_col(self.df_tsne_file)[0]
        self.cell_barcode = [item.split("\t")[0] for item in self.cell_barcode]
        del self.cell_barcode[0]

    @staticmethod
    @utils.add_log
    def get_chunk_df(df_peak, df_fragments):
        index_res = set()
        for ch in set(df_peak.chr):
            df_peak_chr = df_peak[df_peak["chr"] == ch]
            df_fragment_chr = df_fragments[df_fragments["chr"] == ch]
            for _, data_peak in df_peak_chr.iterrows():
                frag_overlap_peak = df_fragment_chr[
                    (df_fragment_chr["start"] >= data_peak["start"])
                    & (df_fragment_chr["end"] <= data_peak["end"])
                ]
                index_res.update(set(frag_overlap_peak.index))
        return index_res

    @utils.add_log
    def count_overlap_peak(self):
        """count fragments overlapping peaks"""
        self.df_fragments.sort_values(["chr", "start", "end"], inplace=True)
        self.df_peaks.sort_values(["chr", "start", "end"], inplace=True)

        peaks_count = self.df_peaks.shape[0]
        chunk_size = peaks_count // self.thread + 1
        df_peak_list = [
            self.df_peaks.iloc[chunk_size * i : chunk_size * (i + 1), :]
            for i in range(self.thread)
        ]
        df_fragment_list = [
            self.df_fragments[self.df_fragments["chr"].isin(set(df_peak.chr))]
            for df_peak in df_peak_list
        ]

        with Pool(self.thread) as p:
            results = p.starmap(Cells.get_chunk_df, zip(df_peak_list, df_fragment_list))

        final_index = set(itertools.chain.from_iterable(results))
        final_df = self.df_fragments[self.df_fragments.index.isin(final_index)]

        final_df_count = final_df.groupby("barcode", as_index=False).agg(
            {"count": "sum"}
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
            name="Estimated Number of Cells",
            value=len(self.cell_barcode),
            help_info="The total number of barcodes identified as cells.",
        )

        raw_reads = self.get_slot_key(
            slot="metrics",
            step_name="barcode",
            key="Raw Reads",
        )
        self.add_metric(
            name="Mean raw read pairs per cell",
            value=int(raw_reads / cell_num),
            help_info="Total number of read pairs divided by the number of cell barcodes.",
        )

        total_fragments = sum(self.df_barcode.fragments)
        cell_fragments = sum(self.df_cell_barcode.fragments)
        self.add_metric(
            name="Fraction of high-quality fragments in cells",
            value=f"{round(cell_fragments / total_fragments * 100, 2)}%",
            help_info="Fraction of high-quality fragments with a valid barcode that are associated with cell-containing partitions.",
        )

        frac_peak = self.df_cell_barcode["frac_peak"].mean()
        self.add_metric(
            name="Fraction of Fragments Overlap with Peaks in Cells",
            value=f"{round(frac_peak * 100, 2)}%",
            help_info="The proportion of fragments in a cell to overlap with a peak.",
        )

        frac_promoter = self.df_cell_barcode["frac_promoter"].mean()
        self.add_metric(
            name="Fraction of fragments overlap with promoter in cells",
            value=f"{round(frac_promoter * 100, 2)}%",
            help_info="The proportion of fragments in a cell to overlap with a promoter sequence.",
        )

        self.add_metric(
            name="Median high-quality fragments per cell",
            value=int(np.median(self.df_cell_barcode.fragments)),
            help_info="The median number of high-quality fragments per cell barcode",
        )

        self.df_barcode.sort_values(by="overlap_peaks", ascending=False, inplace=True)
        self.df_barcode.iloc[0, -1] = True
        self.add_data(
            chart=get_plot_elements.plot_barcode_rank(
                self.df_barcode, log_uniform=False
            )
        )
