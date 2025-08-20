import subprocess
import pandas as pd 
import math
from celescope.tools.plotly_plot import Peak_plot, Tsne_plot, Frag_dis_plot
from celescope.tools import utils
from celescope.tools.step import Step, s_common


def get_opts_analysis(parser, sub_program):

    if sub_program:
        s_common(parser)
        parser.add_argument('--analysis_dir', help='analysis directory', required=True)
        parser.add_argument('--filtered_peak_count', help='filtered peak count h5 file', required=True)
    return parser


class Analysis(Step):
    """Analysis
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        
        self.filtered_peak_count = args.filtered_peak_count
        self.analysis_dir = args.analysis_dir
        self.cell_qc_metrics = f"{self.analysis_dir}/cell_qc_metrics.tsv"
        self.peak_res =  f"{self.analysis_dir}/peak/{self.sample}_final_peaks.bed"
        self.df_tsne_file = f"{self.analysis_dir}/tsne_coord.tsv"
        self.fragment = f"{self.analysis_dir}/fragments_corrected_dedup_count.tsv.gz*"
        self.out = f"{self.outdir}/../outs"

    @staticmethod
    def cell_label(df):
        if df["cell_called"] == True:
            return "Cells"
        else:
            return "Non-cells"
    
    @utils.add_log
    def add_metrics(self):
        """plot and add metrics.
        """
        df = pd.read_csv(self.cell_qc_metrics, sep='\t')
        df['cell_called'] = df.apply(self.cell_label, axis=1)
        
        Peakplot = Peak_plot(df=df).get_plotly_div()
        self.add_data(Peak_plot=Peakplot)

        Fragdisplot = Frag_dis_plot(df=df).get_plotly_div()
        self.add_data(Frag_dis_plot=Fragdisplot)

        df_peak = pd.read_csv(self.peak_res, sep='\t', header=None)
        total_peak = df_peak.shape[0]
        self.add_metric(
            name = 'Number of peaks',
            value = format(int(total_peak), ','),
            help_info = 'Total number of peaks on primary contigs either detected by the pipeline or input by the user.'           
        )
        
        df_tsne = pd.read_csv(self.df_tsne_file, sep='\t')
        df_tsne = df_tsne.rename(columns={"Unnamed: 0": "barcode"})
        df_tsne = pd.merge(df_tsne, df)
        tsne_cluster = Tsne_plot(df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        df_tsne['log10 Fragments'] = df_tsne['fragments'].apply(lambda x: math.log10(x))
        tsne_fragment = Tsne_plot(df_tsne, 'log10 Fragments', discrete=False).get_plotly_div()
        self.add_data(tsne_fragment=tsne_fragment)
    
    def cp_files(self):
        """copy files"""
        files = [self.df_tsne_file, self.peak_res, self.cell_qc_metrics, f"{self.analysis_dir}/peak/*.h5", self.fragment]
        for file in files:
            cmd = f"cp {file} {self.outdir}"
            subprocess.check_call(cmd, shell=True)
        
        utils.check_mkdir(self.out)
        files = [self.cell_qc_metrics, self.filtered_peak_count, self.fragment]
        for file in files:
            cmd = f"cp {file} {self.out}"
            subprocess.check_call(cmd, shell=True)
            
    def run(self):
        self.add_metrics()
        self.cp_files()


def analysis(args):
    with Analysis(args) as runner:
        runner.run()