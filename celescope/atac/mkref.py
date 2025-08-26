import subprocess
import os
import configparser

import celescope
from celescope.tools import utils
from celescope.__init__ import __VERSION__

GENOME_CONFIG = "celescope_atac_genome.config"
TOOLS_DIR = os.path.dirname(celescope.tools.__file__)


class Mkref:
    """
    ## Features

    - Make atac reference.

    ## Usage
    ```
    celescope atac mkref --fasta genome.fa --gtf genes.gtf
    ```

    """

    def __init__(self, args):
        self.fasta = args.fasta
        self.gtf = args.gtf
        self.outdir = f"{os.getcwd()}"
        self.files = {}
        self.meta = {}
        self.files["fasta"] = self.fasta
        self.files["gtf"] = self.gtf
        self.meta["genome_type"] = "atac"
        self.meta["celescope_version"] = __VERSION__

    def __call__(self):
        self.rename()
        self.build_index()
        self.gtf_to_genepred()
        self.gen_promoter()
        self.write_config()

    @utils.add_log
    def rename(self):
        """
        rename to genome.fa and gene.gtf.
        """
        cmd_fa = f"mv {self.fasta} genome.fa"
        cmd_gtf = f"mv {self.gtf} gene.gtf"
        for cmd in [cmd_fa, cmd_gtf]:
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def build_index(self):
        cmd = "chromap -i -r genome.fa -o genome.index"
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def gtf_to_genepred(self):
        cmd = (
            f"{TOOLS_DIR}/gtfToGenePred -genePredExt -geneNameAsName2 gene.gtf gene.tmp"
        )
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def gen_promoter(self):
        cmd = """awk '{if($4>=2000) print $2"\t"$4-2000"\t"$4+2000"\t"$1"\t"$12"\t"$3}' gene.tmp > promoter.bed"""
        subprocess.check_call(cmd, shell=True)

    def write_config(self):
        config = configparser.ConfigParser()
        config.optionxform = str
        config["files"] = self.files
        config["meta"] = self.meta

        with open(GENOME_CONFIG, "w") as config_handle:
            config.write(config_handle)


def mkref(args):
    runner = Mkref(args)
    runner()


def get_opts_mkref(parser, sub_program=True):
    if sub_program:
        parser.add_argument("--fasta", help="Required. fasta file name.", required=True)
        parser.add_argument("--gtf", help="Required. Gtf file name.", required=True)
    return parser
