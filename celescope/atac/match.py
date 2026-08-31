import json
import pandas as pd
import numpy as np
from celescope.tools.matrix import CountMatrix
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools import utils
from celescope.tools.analysis_wrapper import Analysis as Tools_analysis
from celescope.tools.step import Step, s_common


def get_opts_match(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument(
            "--match_dir", help="Matched scRNA-seq directory", required=True
        )
        parser.add_argument(
            "--matrix_file", help="Matrix directory path.", required=True
        )
    return parser


class Cells_metrics(Step):
    @utils.add_log
    def add_cells_metrics(
        self,
        n_cells,
        fraction_reads_in_cells,
        mean_used_reads_per_cell,
        median_umi_per_cell,
        total_genes,
        median_genes_per_cell,
        saturation,
        valid_reads,
    ):
        self.add_metric(
            name="Estimated Number of Cells",
            value=n_cells,
            help_info="The number of barcodes considered as cell-associated.",
        )

        self.add_metric(
            name="Fraction Reads in Cells",
            value=f"{round(fraction_reads_in_cells * 100, 2)}%",
            help_info="the fraction of uniquely-mapped-to-transcriptome reads with cell-associated barcodes",
        )

        mean_reads_per_cell = valid_reads // n_cells
        self.add_metric(
            name="Mean Reads per Cell",
            value=mean_reads_per_cell,
            help_info="the number of Valid Reads divided by Estimated Number of Cells",
        )

        self.add_metric(
            name="Mean Used Reads per Cell",
            value=mean_used_reads_per_cell,
            help_info="the number of uniquely-mapped-to-transcriptome reads per cell-associated barcode",
        )

        self.add_metric(
            name="Median UMI per Cell",
            value=median_umi_per_cell,
            help_info="the median number of UMI counts per cell-associated barcode",
        )

        self.add_metric(
            name="Total Genes",
            value=total_genes,
            help_info="the number of genes with at least one UMI count in any cell",
        )

        self.add_metric(
            name="Median Genes per Cell",
            value=median_genes_per_cell,
            help_info="the median number of genes detected per cell-associated barcode",
        )

        self.add_metric(
            name="Saturation",
            value=f"{round(saturation * 100, 2)}%",
            help_info="the fraction of read originating from an already-observed UMI. ",
        )

    def run(self):
        pass


class Cells_rna(Cells_metrics):
    """Run rna cells and rna analysis to keep the cell numbers of RNA and ATAC consistent"""

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # in
        self.match_dir = args.match_dir
        self.rna_json = f"{self.match_dir}/.metrics.json"
        self.filtered_matrix = args.matrix_file
        self.counts_file = f"{self.match_dir}/outs/counts.tsv"

        # out
        self.filter_counts_files = f"{self.outdir}/counts.tsv"

    @utils.add_log
    def rna_mapping_metrics(self):
        with open(self.rna_json, "r", encoding="utf-8") as file:
            data = json.load(file)
        mapping_metrics = data["mapping_summary"]

        self.add_metric(
            name="genome",
            value=mapping_metrics["Genome"],
        )

        name = "Reads Mapped To Unique Loci"
        self.add_metric(
            name=name,
            value=f"{mapping_metrics[name]}%",
            help_info="Reads that mapped uniquely to the genome.",
        )
        name = "Reads Mapped To Multiple Loci"
        self.add_metric(
            name=name,
            value=f"{mapping_metrics[name]}%",
            help_info="Reads that mapped to multiple loci in the genome",
        )
        name = "Reads Mapped Uniquely To Transcriptome"
        self.add_metric(
            name=name,
            value=f"{mapping_metrics[name]}%",
            help_info="Reads that mapped to a unique gene in the transcriptome. These reads are used for UMI counting.",
        )
        name = "Mapped Reads Assigned To Exonic Regions"
        self.add_metric(
            name=name,
            value=f"{mapping_metrics[name]}%",
            help_info="Reads that assigned to exonic regions of genes",
        )
        name = "Mapped Reads Assigned To Intronic Regions"
        self.add_metric(
            name=name,
            value=f"{mapping_metrics[name]}%",
            help_info="Reads that assigned to intronic regions of genes",
        )
        name = "Mapped Reads Assigned To Intergenic Regions"
        self.add_metric(
            name=name,
            value=f"{mapping_metrics[name]}%",
            help_info="Reads that can not be assigned to a gene will be considered as intergenic reads.",
        )
        name = "Mapped Reads Assigned Antisense To Gene"
        self.add_metric(
            name=name,
            value=f"{mapping_metrics[name]}%",
            help_info="Reads that assigned to the opposite strand of genes",
        )

    @utils.add_log
    def rna_cell_metrics(self, filtered):
        with open(self.rna_json, "r", encoding="utf-8") as file:
            data = json.load(file)
        cells_metrics = data["cells_summary"]

        bcs = filtered.get_barcodes()
        n_cells = len(bcs)

        df_counts = pd.read_csv(self.counts_file, index_col=0, header=0, sep="\t")
        df_counts.index = df_counts.index.str.replace("_", "")
        reads_total = df_counts["countedU"].sum()
        reads_cell = df_counts.loc[bcs, "countedU"].sum()
        fraction_reads_in_cells = float(reads_cell / reads_total)
        mean_used_reads_per_cell = int(reads_cell // len(bcs))
        median_umi_per_cell = int(df_counts.loc[bcs, "UMI"].median())

        bc_geneNum, total_genes = filtered.get_bc_geneNum()
        median_genes_per_cell = int(np.median(list(bc_geneNum.values())))

        saturation = cells_metrics["Saturation"] / 100
        df_counts.loc[:, "mark"] = "UB"
        df_counts.loc[bcs, "mark"] = "CB"
        df_counts.to_csv(self.filter_counts_files, sep="\t", index=True)
        valid_reads = (
            cells_metrics["Mean Reads per Cell"]
            * cells_metrics["Estimated Number of Cells"]
        )
        self.add_cells_metrics(
            n_cells,
            fraction_reads_in_cells,
            mean_used_reads_per_cell,
            median_umi_per_cell,
            total_genes,
            median_genes_per_cell,
            saturation,
            valid_reads,
        )
        self.add_data(
            chart=get_plot_elements.plot_barcode_rank(self.filter_counts_files)
        )

    def run(self):
        if self.match_dir != "None":
            self.rna_mapping_metrics()
            filtered = CountMatrix.from_matrix_dir(self.filtered_matrix)
            self.rna_cell_metrics(filtered)


class Analysis_rna(Tools_analysis):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)


def match(args):
    with Cells_rna(args) as runner:
        runner.run()

    with Analysis_rna(args) as runner:
        runner.run()
