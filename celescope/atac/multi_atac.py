from celescope.tools.multi import Multi, TOOLS_DIR
from celescope.atac.__init__ import __ASSAY__
from celescope.atac.atac import __SUB_STEPS__


class Multi_atac(Multi):
    """
    ## Usage
    ```
    multi_atac \\
        --mapfile mapfile \\
        --thread 8 \\
        --chemistry atac \\
        --reference /path/references \\
        --giggleannotation /path/annotations/giggle.all \\
        --mod shell
    ``` 
    """
    def atac(self, sample):
        step = "atac"
        cmd_line = self.get_cmd_line(step, sample)
        input_path = f'{self.outdir_dic[sample]["barcode"]}'
        cmd = (
            f'{cmd_line} '
            f'--input_path {input_path} '
        )
        self.process_cmd(cmd, step, sample, m=30, x=self.args.thread)

    def analysis(self, sample):
        step = 'analysis'
        feature_matrix = f'{self.outdir_dic[sample]["atac"]}/unfiltered_feature_matrix.rds'
        fragments = f'{self.outdir_dic[sample]["atac"]}/fragments.bed'
        tss_data = f'{self.outdir_dic[sample]["atac"]}/scPipe_atac_stats/tss_plot_data.csv'
        cell_qc_metrics = f'{self.outdir_dic[sample]["atac"]}/cell_qc_metrics.csv'
        binary_matrix = f'{self.outdir_dic[sample]["atac"]}/binary_matrix.rds'
        sce_rds = f'{self.outdir_dic[sample]["atac"]}/scPipe_atac_SCEobject.rds'
        peak_res = f'{self.outdir_dic[sample]["atac"]}/NA_peaks.narrowPeak'
        
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--feature_matrix {feature_matrix} '
            f'--fragments {fragments} '
            f'--tss_data {tss_data} '
            f'--cell_qc_metrics {cell_qc_metrics} '
            f'--binary_matrix {binary_matrix} '
            f'--sce_rds {sce_rds} '
            f'--peak_res {peak_res} '
        )
        self.process_cmd(cmd, step, sample, m=10, x=1)

    def merge_report(self):
        step = "merge_report"
        _index = self.STEPS.index('atac') + 1
        steps_str = ",".join(self.STEPS[:_index] + __SUB_STEPS__ + self.STEPS[_index:])
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
    multi = Multi_atac(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()