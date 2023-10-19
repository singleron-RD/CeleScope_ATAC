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
    parser.add_argument('--reference ', help='Genome reference fasta file', required=True)
    parser.add_argument('--giggleannotation ', help='Path of the giggle annotation file', required=True)
    parser.add_argument('--species ', choices=['GRCh38', 'GRCm38'] help='GRCh38 for human, GRCm38 for mouse', required=True)

    if sub_program:
        s_common(parser)
        parser.add_argument('--input_path', help='input_path from Barcode step.', required=True)
    return parser


class ATAC(Step):
    """
    ## Features  
    - Run ATAC.
    ## Output
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.giggleannotation = args.giggleannotation
        self.input_path = args.input_path
        self.species = args.species
        self.reference = args.reference
        self.outdir = os.path.abspath(self.outdir)
    
    @utils.add_log
    def run_maestro(self):
        """run maestro
        """
        # Step 1. Configure the MAESTRO workflow
        cmd = (
            f"MAESTRO scatac-init --input_path {self.input_path} "
            f"--gzip --species {self.species} --platform 10x-genomics --format fastq --mapping chromap "
            f"--giggleannotation {self.giggleannotation} "
            f"----fasta {self.reference}/{self.species}_genome.fa "
            f"--index {self.reference}/{self.species}_chromap.index "
            f"--whitelist {self.whitelist}"
            f"--cores {self.thread} --directory {self.outdir} "
            f"--annotation --method RP-based --signature human.immune.CIBERSORT "
            f"--rpmodel Enhanced "
            f"--peak_cutoff 100 --count_cutoff 1000 --frip_cutoff 0.2 --cell_cutoff 50 "
        )
        subprocess.check_call(cmd, shell=True)
        
        # Step 2. Configure samples.json file
        cmd = (
            f"MAESTRO samples-init --assay_type scatac --platform 10x-genomics --data_type fastq --data_dir {self.input_path}"
        )
        
        cwd = os.getcwd()
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
        # change dir back 
        os.chdir(cwd)
        
        # Step 3. Run snakemake pipeline
        cmd = (
            f"snakemake -j {self.thread}"
        )
        subprocess.check_call(cmd, shell=True)
    
    
    def run(self):
        self.run_maestro()


def atac(args):
    with ATAC(args) as runner:
        runner.run()
    
    # with Mapping(args) as runner:
    #     runner.run()
        
    # with Cells(args) as runner:
    #     runner.run()
        
        
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
         
        