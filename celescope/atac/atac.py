import subprocess
import pandas as pd 
import numpy as np 
import os
from celescope.__init__ import ROOT_PATH
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.tools import get_plot_elements
from celescope.tools.plotly_plot import Insert_plot


__SUB_STEPS__ = ['mapping', 'cells']


def get_opts_atac(parser, sub_program):
    parser.add_argument('--reference', help='Genome reference fasta file', required=True)
    parser.add_argument('--giggleannotation', help='Path of the giggle annotation file', required=True)
    parser.add_argument('--species', choices=['GRCh38', 'GRCm38'], help='GRCh38 for human, GRCm38 for mouse', required=True)
    parser.add_argument(
        '--signature',
        help='Cell signature file used to annotate cell types',
        default='human.immune.CIBERSORT',
        choices=['human.immune.CIBERSORT', 'mouse.brain.ALLEN', 'mouse.all.facs.TabulaMuris', 'mouse.all.droplet.TabulaMuris'],
    )
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
        self.input_path = os.path.abspath(args.input_path)
        self.species = args.species
        self.reference = args.reference
        self.outdir = os.path.abspath(self.outdir)
        self.signature = args.signature
        self.whitelist = f"{ROOT_PATH}/data/chemistry/atac/857K-2023.txt"
    
    @utils.add_log
    def run_maestro(self):
        """run maestro
        """
        # Step 1. Configure the MAESTRO workflow
        cmd = (
            f"MAESTRO scatac-init --input_path {self.input_path} "
            f"--species {self.species} --platform 10x-genomics --format fastq --mapping chromap "
            f"--giggleannotation {self.giggleannotation} "
            f"--fasta {self.reference}/{self.species}_genome.fa "
            f"--index {self.reference}/{self.species}_chromap.index "
            f"--cores {self.thread} --directory {self.outdir} "
            f"--annotation --method RP-based --signature {self.signature} "
            f"--rpmodel Enhanced "
            f"--peak_cutoff 100 --count_cutoff 1000 --frip_cutoff 0.2 --cell_cutoff 50 "
            f"--whitelist {self.whitelist} "
        )
        subprocess.check_call(cmd, shell=True)
        
        # Step 2. Configure samples.json file
        cmd = (
            f"MAESTRO samples-init --assay_type scatac --platform 10x-genomics --data_type fastq --data_dir {self.input_path}"
        )
        
        cwd = os.getcwd()
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
        
        # Step 3. Run snakemake pipeline
        cmd = (
            f"snakemake -j {self.thread}"
        )
        subprocess.check_call(cmd, shell=True)
        
        # change dir back 
        os.chdir(cwd)
    
    
    def run(self):
        self.run_maestro()


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
        
        self.cell_barcode = utils.read_one_col(f"{self.outdir}/Result/QC/{self.sample}/{self.sample}_scATAC_validcells.txt")[0]
        self.df_barcode = pd.read_csv(f"{self.outdir}/Result/QC/{self.sample}/singlecell.txt",
                                      header=None, sep='\t', names=["barcode", "fragments", "fragments_overlapping_promoter"])
        self.df_cell_barcode = self.df_barcode[self.df_barcode["barcode"].isin(self.cell_barcode)]
        self.df_fragments = pd.read_csv(f"{self.outdir}/Result/Mapping/{self.sample}/fragments_corrected_count_sortedbybarcode.tsv",
                                      header=None, sep='\t', names=["chrom", "chromStart", "chromEnd", "barcode", "count"])
        self.df_peaks = pd.read_csv(f"{self.outdir}/Result/Analysis/{self.sample}/{self.sample}_final_peaks.bed",
                                    header=None, sep='\t', names=["chr", "start", "end"])
        self.df_cell_metrics = f"{self.outdir}/Result/Analysis/{self.sample}/cell_qc_metrics.tsv"
    
    @utils.add_log
    def count_overlap_peak(self):
        """count fragments overlapping peaks
        """
        final_df = pd.DataFrame()
        for _, data_peak in self.df_peaks.iterrows():
            peak_info = dict(data_peak)
            frag_chr = self.df_fragments[self.df_fragments["chr"] == peak_info["chr"]]
            frag_overlap_peak = frag_chr[ (frag_chr["start"] >= peak_info["start"]) & (frag_chr["end"] <= peak_info["end"])] 
            final_df = pd.concat([final_df, frag_overlap_peak])
            final_df_count = final_df.groupby('barcode', as_index=False).agg({'count': 'sum'})
            
            self.df_barcode = pd.merge(self.df_barcode, final_df_count, on='barcode', how='outer').fillna(0)
            self.df_barcode = self.df_barcode.rename(columns={'count': 'overlap_peaks'})
            self.df_barcode['frac_peak'] = round(self.df_barcode['overlap_peaks'] / self.df_barcode['fragments'], 4)
            self.df_barcode['cell_called'] = self.df_barcode['barcode'].apply(lambda x: True if x in self.cell_barcode else False)
            self.df_barcode.to_csv(self.df_cell_metrics, sep='\t', index=False)
            
    def run(self):
        self.count_overlap_peak()
        #self.df_barcode = pd.read_csv(self.df_cell_metrics, sep='\t')
        
        cell_num = len(self.cell_barcode)
        self.add_metric(
            name = 'Estimated Number of Cells',
            value = len(self.cell_barcode),
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

        total_fragments = sum(self.df_barcode.fragments)
        cell_fragments = sum(self.df_cell_barcode.fragments)
        self.add_metric(
            name = 'Fraction of high-quality fragments in cells',
            value = f'{round(cell_fragments / total_fragments * 100, 2)}%',
            help_info = 'Fraction of high-quality fragments with a valid barcode that are associated with cell-containing partitions.'           
        )

        frac_peak = self.df_barcode[self.df_barcode['cell_called']==True]['frac_peak'].mean()
        self.add_metric(
            name = 'Fraction of Fragments Overlap with Peaks in Cells',
            value = f'{round(frac_peak * 100, 2)}%',
            help_info = 'The proportion of fragments in a cell to overlap with a peak.'           
        )
        
        total_promoter = sum(self.df_barcode.fragments_overlapping_promoter)
        cell_promoter = sum(self.df_cell_barcode.fragments_overlapping_promoter)
        self.add_metric(
            name = 'Fraction of fragments overlap with promoter in cells',
            value = f'{round(cell_promoter / total_promoter * 100, 2)}%',
            help_info = 'The proportion of fragments in a cell to overlap with a promoter sequence.'           
        )

        self.add_metric(
            name = 'Median high-quality fragments per cell',
            value = np.median(self.df_cell_barcode.fragments),
            help_info = 'The median number of high-quality fragments per cell barcode'           
        )
        
        self.df_barcode.sort_values(by='overlap_peaks', ascending=False, inplace=True)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.df_barcode, log_uniform=False))
        
        
class Mapping(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.df_mapping = pd.read_csv(f"{self.outdir}/Result/Mapping/{self.sample}/fragments_corrected_count_sortedbybarcode.tsv",
                                      header=None, sep='\t', names=["chrom", "chromStart", "chromEnd", "barcode", "count"])
    
    def run(self):

        valid_reads = self.get_slot_key(
            slot='metrics',
            step_name='barcode',
            key='Valid Reads',
        )
        
        self.df_mapping["size"] = self.df_mapping["chromEnd"] - self.df_mapping["chromStart"]
        
        self.add_metric(
            name = 'Confidently Unique mapped read pairs',
            value = sum(self.df_mapping['count']),
            total = valid_reads,
            help_info = 'Fraction of Unique mapped read pairs'
        )

        self.add_metric(
            name = 'Fragments in nucleosome-free regions',
            value = sum(self.df_mapping[self.df_mapping['size']<=124]['count']),
            total = sum(self.df_mapping['count']),
            help_info = 'Fraction of high-quality fragments smaller than 124 basepairs.'
        )

        self.add_metric(
            name = 'Fragments flanking a single nucleosome',
            value = sum(self.df_mapping[ (self.df_mapping['size']>124) & (self.df_mapping['size']<=296)]['count']),
            total = sum(self.df_mapping['count']),
            help_info = 'Fraction of high-quality fragments between 124 and 296 basepairs.'
        )
        self.df_mapping = self.df_mapping[self.df_mapping['size'] < 800]
        Insertplot = Insert_plot(df=self.df_mapping).get_plotly_div()
        self.add_data(Insert_plot=Insertplot)
         
        