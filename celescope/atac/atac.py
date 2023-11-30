import subprocess
import pandas as pd 
import numpy as np 
import os
import itertools
from multiprocessing import Pool
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
    parser.add_argument('--peak_cutoff', type=int, help='Minimum number of peaks included in each cell', default=100)
    parser.add_argument('--count_cutoff', type=int, help='Cutoff for the number of count in each cell', default=1000)
    parser.add_argument('--frip_cutoff', type=float, help='Cutoff for fraction of reads in promoter in each cell', default=0.2)
    parser.add_argument('--cell_cutoff', type=int,  help='Minimum number of cells covered by each peak', default=10)
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
        
        # cut-off
        self.peak_cutoff = args.peak_cutoff
        self.count_cutoff = args.count_cutoff
        self.frip_cutoff = args.frip_cutoff
        self.cell_cutoff = args.cell_cutoff
    
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
            f"--peak_cutoff {self.peak_cutoff} --count_cutoff {self.count_cutoff} --frip_cutoff {self.frip_cutoff} --cell_cutoff {self.cell_cutoff} "
            f"--whitelist {self.whitelist} "
            f"--shortpeak "
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


class Cells(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        
        self.df_barcode = pd.read_csv(f"{self.outdir}/Result/QC/{self.sample}/singlecell.txt",
                                      header=None, sep='\t', names=["barcode", "fragments", "fragments_overlapping_promoter"])
        self.df_fragments = pd.read_csv(f"{self.outdir}/Result/Mapping/{self.sample}/fragments_corrected_count_sortedbybarcode.tsv",
                                      header=None, sep='\t', names=["chr", "start", "end", "barcode", "count"])
        self.df_peaks = pd.read_csv(f"{self.outdir}/Result/Analysis/{self.sample}/{self.sample}_final_peaks.bed",
                                    header=None, sep='\t', names=["chr", "start", "end"])
        
        # out
        self.sce_rds = f"{self.outdir}/Result/Analysis/{self.sample}/{self.sample}_scATAC_Object.rds"
        self.df_cell_metrics = f"{self.outdir}/Result/Analysis/{self.sample}/cell_qc_metrics.tsv"
        self.meta_data = f"{self.outdir}/Result/Analysis/{self.sample}/meta.csv"
    
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

    @staticmethod
    @utils.add_log
    def get_chunk_df(df_peak, df_fragments):
        index_res = set()
        for ch in set(df_peak.chr):
            df_peak_chr = df_peak[df_peak['chr'] == ch]
            df_fragment_chr = df_fragments[df_fragments['chr']== ch]
            for _, data_peak in df_peak_chr.iterrows():
                frag_overlap_peak = df_fragment_chr[ (df_fragment_chr["start"] >= data_peak["start"]) & (df_fragment_chr["end"] <= data_peak["end"])]
                index_res.update(set(frag_overlap_peak.index))
        return index_res
    
    @utils.add_log
    def count_overlap_peak(self):
        """count fragments overlapping peaks
        """
        self.cell_barcode = utils.read_one_col(self.meta_data)[0]
        del self.cell_barcode[0]
        
        self.df_fragments.sort_values(['chr','start','end'], inplace=True)
        self.df_peaks.sort_values(['chr','start','end'], inplace=True)
        
        peaks_count = self.df_peaks.shape[0]
        chunk_size = peaks_count // self.thread + 1
        df_peak_list = [self.df_peaks.iloc[ chunk_size*i: chunk_size*(i+1), :] for i in range(self.thread)]
        df_fragment_list = [self.df_fragments[self.df_fragments['chr'].isin(set(df_peak.chr))] for df_peak in df_peak_list]
        
        with Pool(self.thread) as p:
            results = p.starmap(Cells.get_chunk_df, zip(df_peak_list, df_fragment_list))
        
        final_index = set(itertools.chain.from_iterable(results))
        final_df = self.df_fragments[self.df_fragments.index.isin(final_index)]
        
        final_df_count = final_df.groupby('barcode', as_index=False).agg({'count': 'sum'})
        self.df_barcode = pd.merge(self.df_barcode, final_df_count, on='barcode', how='outer').fillna(0)
        self.df_barcode = self.df_barcode.rename(columns={'count': 'overlap_peaks'})
        self.df_barcode['frac_peak'] = round(self.df_barcode['overlap_peaks'] / self.df_barcode['fragments'], 4)
        self.df_barcode['cell_called'] = self.df_barcode['barcode'].apply(lambda x: True if x in self.cell_barcode else False)
        self.df_barcode.to_csv(self.df_cell_metrics, sep='\t', index=False)
            
    def run(self):
        self.gen_plot_data()
        self.count_overlap_peak()
        
        self.df_cell_barcode = self.df_barcode[self.df_barcode["barcode"].isin(self.cell_barcode)]
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
            value = int(np.median(self.df_cell_barcode.fragments)),
            help_info = 'The median number of high-quality fragments per cell barcode'           
        )
        
        self.df_barcode.sort_values(by='overlap_peaks', ascending=False, inplace=True)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.df_barcode, log_uniform=False))
         
        