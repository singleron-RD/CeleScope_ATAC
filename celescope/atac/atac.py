import subprocess
import pandas as pd 
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import ROOT_PATH


def get_opts_atac(parser, sub_program):
    parser.add_argument('--reference', help='reference path', required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--r1', help='R1 reads from step Barcode.', required=True)
        parser.add_argument('--r2', help='R2 reads from step Barcode.', required=True)
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
        self.reference = args.reference
    
    @utils.add_log
    def run_scpipe(self):
        """run scpipe
        """
        cmd = (
            f"Rscript {ROOT_PATH}/atac/scpipe.R "
            f"--sample {self.sample} "
            f"--r1 {self.r1} "
            f"--r2 {self.r2} "
            f"--thread {self.thread} "
            f"--outdir {self.outdir} "
            f"--reference {self.reference}"
        )
        subprocess.check_call(cmd, shell=True)
    
    def run(self):
        self.run_scpipe()


def atac(args):
    
    with ATAC(args) as runner:
        runner.run()