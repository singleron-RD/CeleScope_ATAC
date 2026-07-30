import subprocess
import re
import os
import celescope.tools
from celescope.tools import utils
from celescope.tools.step import Step, s_common

TOOLS_DIR = os.path.dirname(celescope.tools.__file__) + "/multi_arc"


def get_opts_cells(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument(
            "--match_dir", help="Matched scRNA-seq directory", required=True
        )
        parser.add_argument(
            "--analysis_dir", help="atac analysis directory", required=True
        )
    return parser


class Cells(Step):
    """Run rna cells and rna analysis to keep the cell numbers of RNA and ATAC consistent"""

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # in
        self.match_dir = args.match_dir
        self.df_tsne_file = f"{args.analysis_dir}/tsne_coord.tsv"
        self.cmd_line = f"{self.outdir}/{self.sample}_cmd_line"

        # out
        self.atac_barcode = f"{self.outdir}/barcodes.tsv"

    @utils.add_log
    def write_barcode(self):
        cell_barcode = utils.read_one_col(self.df_tsne_file)[0]
        cell_barcode = [item.split("\t")[0] for item in cell_barcode]
        del cell_barcode[0]
        cell_barcode = [
            item[:9] + "_" + item[9:18] + "_" + item[18:] for item in cell_barcode
        ]
        with open(self.atac_barcode, "w") as fh:
            for barcode in cell_barcode:
                fh.write(f"{barcode}\n")

    @staticmethod
    def get_rna_genomedir(match_dir):
        sjm_file = f"{match_dir}/../sjm/sjm.job"
        pattern = r"--genomeDir\s+([^\s]+)"
        with open(sjm_file, "r") as f:
            for line in f:
                match = re.search(pattern, line)
                if match:
                    genome_dir = match.group(1)
                    break

        return genome_dir

    @utils.add_log
    def rna_cells(self):
        rna_sample = self.match_dir.split("/")[-1]
        rna_outdir = self.match_dir
        genome_dir = Cells.get_rna_genomedir(self.match_dir)
        cmd = (
            "source activate celescope3.0.0\n"
            "celescope rna cells "
            f"--barcode {os.path.abspath(self.atac_barcode)} "
            f"--outdir {rna_outdir}/cells "
            f"--sample {rna_sample}\n"
            "celescope rna analysis "
            f"--outdir {rna_outdir}/02.analysis "
            f"--sample {rna_sample} "
            f"--genomeDir {genome_dir} "
            f"--matrix_file {rna_outdir}/outs/filtered "
        )
        with open(self.cmd_line, "w") as f:
            f.write(cmd)

        cwd = os.getcwd()
        # avoid celescope.tools.matrix.from_matrix_dir error
        os.chdir(f"{rna_outdir}/..")
        subprocess.check_call(cmd, shell=True)
        # change dir back to avoid error
        os.chdir(cwd)

    @utils.add_log
    def merge_report(self):
        cmd = f"python {TOOLS_DIR}/merge.py --rna {self.match_dir} --atac {self.outdir}/.. --outdir {self.outdir}/.."
        subprocess.check_call(cmd, shell=True)

    def run(self):
        if self.match_dir != "None":
            self.write_barcode()
            self.rna_cells()
            self.merge_report()


def cells(args):
    cells_obj = Cells(args)
    cells_obj.run()
