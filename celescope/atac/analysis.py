import subprocess
import pandas as pd 
from celescope.tools.plotly_plot import Insert_plot, Tss_plot, Peak_plot
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import ROOT_PATH


def get_opts_analysis(parser, sub_program):

    if sub_program:
        s_common(parser)
        parser.add_argument('--feature_matrix', help='unfiltered feature matrix rds file.', required=True)
        parser.add_argument('--fragments', help='fragments bed file.', required=True)
        parser.add_argument('--tss_data', help='tss plot data csv file.', required=True)
        parser.add_argument('--cell_qc_metrics', help='cell qc metrics csv file.', required=True)
        parser.add_argument('--binary_matrix', help='binary matrix rds file.', required=True)
        parser.add_argument('--sce_rds', help='scPipe atac SCEobject rds file.', required=True)
        parser.add_argument('--peak_res', help='peak result excel file.', required=True)
    return parser


class Analysis(Step):
    """Analysis
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        
        self.feature_matrix = args.feature_matrix
        self.fragments = args.fragments
        self.tss_data = args.tss_data
        self.cell_qc_metrics = args.cell_qc_metrics
        self.binary_matrix = args.binary_matrix
        self.sce_rds = args.sce_rds
        self.peak_res = args.peak_res
        
        self.insert_size_data = f'{self.outdir}/insert_size.csv'
        self.cell_umap_info = f'{self.outdir}/cell_umap.csv'

    @utils.add_log
    def gen_plot_data(self):
        """generate data for plot.
        """
        cmd = (
            f"Rscript {ROOT_PATH}/atac/gen_plot_data.R "
            f"--feature_matrix {self.feature_matrix} "
            f"--fragments {self.fragments} "
            f"--insert_size_data {self.insert_size_data} "
            f"--binary_matrix {self.binary_matrix} "
            f"--sce_rds {self.sce_rds} "
            f"--cell_umap_info {self.cell_umap_info} "
        )
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def add_metrics(self):
        """plot and add metrics.
        """
        df_insert = pd.read_csv(self.insert_size_data)
        Insertplot = Insert_plot(df=df_insert).get_plotly_div()
        self.add_data(Insert_plot=Insertplot)
        
        df_tss = pd.read_csv(self.tss_data)
        Tssplot = Tss_plot(df=df_tss).get_plotly_div()
        self.add_data(Tss_plot=Tssplot)
        
        df_cell = pd.read_csv(self.cell_qc_metrics)
        Peakplot = Peak_plot(df=df_cell).get_plotly_div()
        self.add_data(Peak_plot=Peakplot)

        with open(self.peak_res) as fh:
            for line in fh.readlines():
                if '# total tags in treatment' in line:
                    total_peak = line.split(':')[-1].strip()
                    break
                
        self.add_metric(
            name = 'Number of peaks',
            value = format(int(total_peak), ','),
            help_info = 'Total number of peaks on primary contigs either detected by the pipeline or input by the user.'           
        )

        self.add_metric(
            name = 'TSS enrichment score',
            value = round(max(df_tss.agg_tss_scores), 2),
            help_info = 'Maximum value of the transcription-start-site (TSS) profile.The TSS profile is the summed accessibility signal (defined as number of cut sites per base) in a window of 2,000 bases around all the annotated TSSs, normalized by the minimum signal in the window.'           
        )

    def run(self):
        self.gen_plot_data()
        self.add_metrics()


def analysis(args):
    with Analysis(args) as runner:
        runner.run()