import subprocess
import pandas as pd 
import numpy as np 
import os
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import ROOT_PATH
from celescope.tools import get_plot_elements


__SUB_STEPS__ = ['mapping', 'cells']


def get_opts_atac(parser, sub_program):
    parser.add_argument('--reference', help='reference path', required=True)
    parser.add_argument('--organism', help='hg38 for human, mm10 for mouse', required=True)
    parser.add_argument(
        '--rm_pcr_dup', 
        help="Whether or not to remove PCR duplicate. High-quality read pairs that are deemed to be PCR duplicates.",
        action='store_true')
    if sub_program:
        s_common(parser)
        parser.add_argument('--r1', help='R1 reads from step Barcode.', required=True)
        parser.add_argument('--r2', help='R2 reads from step Barcode.', required=True)
        parser.add_argument('--r3', help='R3 reads from step Barcode.', required=True)
    return parser


class ATAC(Step):
    """
    ## Features  
    - Run ATAC.
    ## Output
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.r1 = args.r1
        self.r2 = args.r2
        self.r3 = args.r3
        self.reference = args.reference
        self.organism = args.organism
        self.rm_pcr_dup = args.rm_pcr_dup
        self.outdir = os.path.abspath(self.outdir)
    
    @utils.add_log
    def run_scpipe(self):
        """run scpipe
        """
        cmd = (
            f"Rscript {ROOT_PATH}/atac/scpipe.R "
            f"--sample {self.sample} "
            f"--r1 {self.r1} "
            f"--r2 {self.r2} "
            f"--r3 {self.r3} "
            f"--thread {self.thread} "
            f"--outdir {self.outdir} "
            f"--reference {self.reference} "
            f"--organism {self.organism} "
            f"--rm_pcr_dup {self.rm_pcr_dup} "
        )
        subprocess.check_call(cmd, shell=True)
    
    @utils.add_log
    def cp_report(self):
        cmd = (f"cp {self.outdir}/scPipe_atac_stats/scPipe_atac_report.html {self.sample}")
        os.system(cmd)
    
    def run(self):
        self.run_scpipe()
        #self.cp_report()


def atac(args):
    with ATAC(args) as runner:
        runner.run()
    
    with Mapping(args) as runner:
        runner.run()
        
    with Cells(args) as runner:
        runner.run()
        
        
class Cells(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.df_cell = pd.read_csv(f"{self.outdir}/cell_qc_metrics.csv", index_col=0, sep=',')

    def run(self):
        cell_num = self.df_cell[self.df_cell.cell_called==True].shape[0]
        self.add_metric(
            name = 'Estimated Number of Cells',
            value = cell_num,
            help_info = 'The total number of barcodes identified as cells.'           
        )

        raw_reads = self.get_slot_key(
            slot='metrics',
            step_name='barcode',
            key='Raw Reads',
        )
        self.add_metric(
            name = 'Mean raw read pairs per cell',
            value = int(raw_reads / cell_num),
            help_info = 'Total number of read pairs divided by the number of cell barcodes.'           
        )

        self.add_metric(
            name = 'Fraction of fragments overlap with peaks in cells',
            value = f'{round(np.mean(self.df_cell.frac_peak)* 100, 2)}%',
            help_info = 'The proportion of fragments in a cell to overlap with a peak.'           
        )

        self.add_metric(
            name = 'Fraction of fragments overlap with tss in cells',
            value = f'{round(np.mean(self.df_cell.frac_tss)* 100, 2)}%',
            help_info = 'The proportion of fragments in a cell to overlap with a tss.'           
        )
        
        self.add_metric(
            name = 'Fraction of fragments overlap with enhancer in cells',
            value = f'{round(np.mean(self.df_cell.frac_enhancer)* 100, 2)}%',
            help_info = 'The proportion of fragments in a cell to overlap with a enhancer sequence.'           
        )
        
        self.add_metric(
            name = 'Fraction of fragments overlap with promoter in cells',
            value = f'{round(np.mean(self.df_cell.frac_promoter)* 100, 2)}%',
            help_info = 'The proportion of fragments in a cell to overlap with a promoter sequence.'           
        )
        
        self.add_metric(
            name = 'Fraction of fragments in cells that are mitochondrial',
            value = f'{round(np.mean(self.df_cell.frac_mito)* 100, 2)}%',
            help_info = 'The proportion of fragments in a cell that are mitochondrial.'           
        )

        self.df_cell.sort_values(by='total_frags', ascending=False, inplace=True)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.df_cell, log_uniform=False))
        
        
class Mapping(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.df_mapping = pd.read_csv(f"{self.outdir}/scPipe_atac_stats/stats_file_align.txt", header=0, names=["name", "value"])
        self.mapping_dict = dict(zip(list(self.df_mapping.name), list(self.df_mapping.value)))
    
    def run(self):

        self.add_metric(
            name = 'Confidently mapped read pairs',
            value = self.mapping_dict["Mapped_fragments"],
            total = self.mapping_dict["Total_fragments"],
            help_info = 'Fraction of mapped read pairs'
        )

        self.add_metric(
            name = 'Unique mapped read pairs',
            value = self.mapping_dict["Uniquely_mapped_fragments"],
            total = self.mapping_dict["Total_fragments"],
            help_info = 'Fraction of unique mapped read pairs'
        )

        self.add_metric(
            name = 'Multi mapped read pairs',
            value = self.mapping_dict["Multi_mapping_fragments"],
            total = self.mapping_dict["Total_fragments"],
            help_info = 'Fraction of multi mapped read pairs'
        )
        
        self.add_metric(
            name = 'Unmapped read pairs',
            value = self.mapping_dict["Unmapped_fragments"],
            total = self.mapping_dict["Total_fragments"],
            help_info = 'Fraction of unmapped read pairs'
        )
         
        