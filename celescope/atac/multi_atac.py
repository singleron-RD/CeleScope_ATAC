from celescope.tools.multi import Multi
from celescope.atac.__init__ import __ASSAY__


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


def main():
    multi = Multi_bulk_vdj(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()