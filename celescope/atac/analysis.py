import subprocess
import pandas as pd 
import math
from celescope.tools.plotly_plot import Peak_plot, Umap_plot
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import ROOT_PATH


def get_opts_analysis(parser, sub_program):

    if sub_program:
        s_common(parser)
        parser.add_argument('--cell_qc_metrics', help='cell qc metrics file.', required=True)
        parser.add_argument('--sce_rds', help='scPipe atac SCEobject rds file.', required=True)
        parser.add_argument('--peak_res', help='narrowPeak file.', required=True)
    return parser


class Analysis(Step):
    """Analysis
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.cell_qc_metrics = args.cell_qc_metrics
        self.sce_rds = args.sce_rds
        self.peak_res = args.peak_res
        self.meta_data = f'{self.outdir}/meta.csv'

    @utils.add_log
    def gen_plot_data(self):
        """generate meta-data file with umap coord for plot.
        """
        cmd = (
            f"Rscript {ROOT_PATH}/atac/gen_plot_data.R "
            f"--sce_rds {self.sce_rds} "
            f"--meta_data {self.meta_data} "
        )
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def add_metrics(self):
        """plot and add metrics.
        """
        df = pd.read_csv(self.cell_qc_metrics, sep='\t')
        Peakplot = Peak_plot(df=df).get_plotly_div()
        self.add_data(Peak_plot=Peakplot)

        df_peak = pd.read_csv(self.peak_res, sep='\t', header=None)
        total_peak = df_peak.shape[0]
        self.add_metric(
            name = 'Number of peaks',
            value = format(int(total_peak), ','),
            help_info = 'Total number of peaks on primary contigs either detected by the pipeline or input by the user.'           
        )
        
        df_meta = pd.read_csv(self.meta_data)
        df_meta = df_meta.rename(columns={"Row.names": "barcode", "seurat_clusters": "cluster"})
        df_meta = pd.merge(df_meta, df)
        umap_cluster = Umap_plot(df_meta, 'cluster').get_plotly_div()
        self.add_data(umap_cluster=umap_cluster)

        df_meta['log10 Fragments'] = df_meta['fragments'].apply(lambda x: math.log10(x))
        umap_fragment = Umap_plot(df_meta, 'log10 Fragments', discrete=False).get_plotly_div()
        self.add_data(umap_fragment=umap_fragment)

    def run(self):
        self.gen_plot_data()
        self.add_metrics()


def analysis(args):
    with Analysis(args) as runner:
        runner.run()