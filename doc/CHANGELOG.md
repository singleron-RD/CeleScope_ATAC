## 2024-04-09
 - Add scatac peakcount and scatac filter script to generate h5 matrix file.
 - Replace separator to (":", "-") in h5 file to avoid error of downstream analysis.
 - No longer use Maestro package.

## 2024-02-28
 - Fix `Fraction of fragments overlap with promoter in cells` Metric.
 - Add gtfToGenePred script to generate promoter file.
 - Support customized species analysis.

 ## 2023-12-01
 - Add Saturation Metric `Percent duplicates` .
 - Avoid load ensdb error: Replace Maestro with Seurat to run dim reduction.
 - Modify outs dir.
 
 ## 2023-10-20
 - Do not use giggle to annotate cell type.
 - Use multiprocessing to count overlap peaks.
 - Improve speed and Reduce memory usage.
