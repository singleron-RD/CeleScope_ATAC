import subprocess
import snapatac2 as snap
import pandas as pd
import math
import os
from celescope.tools.plotly_plot import Peak_plot, Tsne_plot, Frag_dis_plot
from celescope.tools import utils
from celescope.tools.step import Step, s_common


def get_opts_analysis(parser, sub_program):
    parser.add_argument(
        "--reference",
        help="Genome reference directory includes fasta, index, promoter bed file.",
        required=True,
    )
    if sub_program:
        s_common(parser)
        parser.add_argument("--analysis_dir", help="analysis directory", required=True)
        parser.add_argument(
            "--filtered_peak_count", help="filtered peak count h5 file", required=True
        )
    return parser


class Analysis(Step):
    """Analysis"""

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.reference = os.path.abspath(args.reference)
        self.filtered_peak_count = args.filtered_peak_count
        self.analysis_dir = args.analysis_dir
        self.fragment = f"{self.analysis_dir}/fragments_corrected_dedup_count.tsv"
        self.cell_qc_metrics = f"{self.analysis_dir}/cell_qc_metrics.tsv"
        self.peak_res = f"{self.analysis_dir}/peak/{self.sample}_final_peaks.bed"
        self.df_tsne_file = f"{self.analysis_dir}/tsne_coord.tsv"
        self.fragment_files = (
            f"{self.analysis_dir}/fragments_corrected_dedup_count.tsv.gz*"
        )
        self.out = f"{self.outdir}/../outs"

    @staticmethod
    def cell_label(df):
        if df["cell_called"] == True:
            return "Cells"
        else:
            return "Non-cells"

    @utils.add_log
    def run_snapatac2(self):
        """run snapatac2"""

        genome = snap.genome.Genome(
            fasta=f"{self.reference}/genome.fa", annotation=f"{self.reference}/gene.gtf"
        )
        data = snap.pp.import_data(
            self.fragment,
            chrom_sizes=genome,
            min_num_fragments=0,
            sorted_by_barcode=False,
            chrM=["chrM", "M", "MT", "Mt"],
        )
        cell_barcode = utils.read_one_col(self.df_tsne_file)[0]
        cell_barcode = [item.split("\t")[0] for item in cell_barcode]
        data_cell = data[data.obs.index.isin(cell_barcode)]
        snap.metrics.tsse(data_cell, genome)
        data_cell.obs.tsse = data_cell.obs.tsse.astype(float)
        tsse = round(max(data_cell.obs["tsse"]), 2)
        self.add_metric(
            name="TSS enrichment score",
            value=tsse,
            help_info="Transcription Start Site (TSS) Enrichment Score.",
        )

    @utils.add_log
    def add_metrics(self):
        """plot and add metrics."""
        df = pd.read_csv(self.cell_qc_metrics, sep="\t")
        df["cell_called"] = df.apply(self.cell_label, axis=1)

        Peakplot = Peak_plot(df=df).get_plotly_div()
        self.add_data(Peak_plot=Peakplot)

        Fragdisplot = Frag_dis_plot(df=df).get_plotly_div()
        self.add_data(Frag_dis_plot=Fragdisplot)

        df_peak = pd.read_csv(self.peak_res, sep="\t", header=None)
        total_peak = df_peak.shape[0]
        self.add_metric(
            name="Number of peaks",
            value=format(int(total_peak), ","),
            help_info="Total number of peaks on primary contigs either detected by the pipeline or input by the user.",
        )

        df_tsne = pd.read_csv(self.df_tsne_file, sep="\t")
        df_tsne = df_tsne.rename(columns={"Unnamed: 0": "barcode"})
        df_tsne = pd.merge(df_tsne, df)
        tsne_cluster = Tsne_plot(df_tsne, "cluster").get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        df_tsne["log10 Fragments"] = df_tsne["fragments"].apply(lambda x: math.log10(x))
        tsne_fragment = Tsne_plot(
            df_tsne, "log10 Fragments", discrete=False
        ).get_plotly_div()
        self.add_data(tsne_fragment=tsne_fragment)

    def cp_files(self):
        """copy files"""
        files = [
            self.df_tsne_file,
            self.peak_res,
            self.cell_qc_metrics,
            f"{self.analysis_dir}/peak/*.h5",
            self.fragment_files,
        ]
        for file in files:
            cmd = f"cp {file} {self.outdir}"
            subprocess.check_call(cmd, shell=True)

        utils.check_mkdir(self.out)
        files = [self.cell_qc_metrics, self.filtered_peak_count, self.fragment_files]
        for file in files:
            cmd = f"cp {file} {self.out}"
            subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_snapatac2()
        self.add_metrics()
        self.cp_files()


def analysis(args):
    with Analysis(args) as runner:
        runner.run()
