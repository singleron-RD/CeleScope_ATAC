from celescope.tools.multi import Multi, TOOLS_DIR
from celescope.atac.__init__ import __ASSAY__
from celescope.atac.atac import __SUB_STEPS__


class Multi_bulk_vdj(Multi):
    """
    ## Usage
    ```
    multi_atac \\
        --mapfile mapfile \\
        --thread 8 \\
        --chemistry atac \\
        --reference /path/GRCh38/fasta/genome.fa \\
        --mod shell
    ``` 
    """
    def atac(self, sample):
        step = "atac"
        cmd_line = self.get_cmd_line(step, sample)
        r1 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_1.fq{self.fq_suffix}'
        r2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq{self.fq_suffix}'
        r3 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_3.fq{self.fq_suffix}'
        cmd = (
            f'{cmd_line} '
            f'--r1 {r1} '
            f'--r2 {r2} '
            f'--r3 {r3} '
        )
        self.process_cmd(cmd, step, sample, m=30, x=self.args.thread)


    def merge_report(self):
        step = "merge_report"
        _index = self.STEPS.index('assemble') + 1
        steps_str = ",".join(self.STEPS[:_index] + __SUB_STEPS__ + self.STEPS[_index:-1])
        samples = ','.join(self.fq_dict.keys())
        app = TOOLS_DIR + '/merge_table.py'
        cmd = (
            f'python {app} --samples {samples} '
            f'--steps {steps_str} --outdir {self.args.outdir}'
        )
        if self.args.rm_files:
            cmd += ' --rm_files'
        self.generate_cmd(cmd, step, sample="")
        for sample in self.fq_dict:
            self.sjm_order += f'order {step} after {self.last_step}_{sample}\n'


def main():
    multi = Multi_bulk_vdj(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()