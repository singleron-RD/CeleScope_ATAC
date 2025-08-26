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

    def barcode(self, sample):
        step = "barcode"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        match_dir = f"{self.col4_dict[sample]}"
        cmd = (
            f"{cmd_line} "
            f"--fq1 {arr[0]} --fq2 {arr[1]} --fq3 {arr[2]} "
            f"--match_dir {match_dir} "
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def atac(self, sample):
        step = "atac"
        cmd_line = self.get_cmd_line(step, sample)
        input_path = f'{self.outdir_dic[sample]["barcode"]}'
        match_dir = f"{self.col4_dict[sample]}"
        cmd = f"{cmd_line} " f"--input_path {input_path} " f"--match_dir {match_dir} "
        self.process_cmd(cmd, step, sample, m=30, x=self.args.thread)

    def analysis(self, sample):
        step = "analysis"
        analysis_dir = f'{self.outdir_dic[sample]["atac"]}'
        filtered_peak_count = (
            f'{self.outdir_dic[sample]["atac"]}/peak/{sample}_filtered_peak_count.h5'
        )

        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f"{cmd_line} "
            f"--analysis_dir {analysis_dir} "
            f"--filtered_peak_count {filtered_peak_count} "
        )
        self.process_cmd(cmd, step, sample, m=10, x=1)

    def merge_report(self):
        step = "merge_report"
        _index = self.STEPS.index("atac") + 1
        steps_str = ",".join(self.STEPS[:_index] + __SUB_STEPS__ + self.STEPS[_index:])
        samples = ",".join(self.fq_dict.keys())
        app = TOOLS_DIR + "/merge_table.py"
        cmd = (
            f"python {app} --samples {samples} "
            f"--steps {steps_str} --outdir {self.args.outdir}"
        )
        if self.args.rm_files:
            cmd += " --rm_files"
        self.generate_cmd(cmd, step, sample="")
        for sample in self.fq_dict:
            self.sjm_order += f"order {step} after {self.last_step}_{sample}\n"


def main():
    multi = Multi_atac(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
