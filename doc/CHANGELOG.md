## [1.6.0] - 2025-08-05
 - Bug fix: Negative count value in peak count file.

## [1.5.2] - 2025-07-23
 - Remove barcode-whitelist when mapping.
 - Replace with scRNA barcode in barcode step.

## [1.5.1] - 2025-07-04
 - update atac-V3.1 pattern.

## [1.5.0] - 2025-05-20
 - Match scRNA barcodes if provide 4th column in mapfile.
 - Delete GGTGGGTGGGTG in 884K-2025 barcode list file.

## [1.4.0] - 2025-01-10
 - add atac-V3 pattern.

## [1.3.0] - 2025-01-10
 - Bug fix: Cannot install gtfToGenePred.
 - Fix plotly version==5.24.1.

## [1.2.0] - 2024-12-10
 - Add mkref step.

## [1.1.0] - 2024-11-29
 - Add --coef for auto peak-cutoff.

## 2024-09-27
 - Add Fragment distribution plot in html.

## 2024-09-06
 - Change the default value of --peak_cutoff to auto to improve cell filter.
 - Add parameter --expected_target_cell_num.
 - Fix error in checkna(X): Missing values found in 'X'.

## 2024-08-07
 - Auto detect chemistry from R1-read.
 - Remove bclist parameter.

## 2024-05-22
 - Use chromap-0.2.6 to avoid mapping fragments error.

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
