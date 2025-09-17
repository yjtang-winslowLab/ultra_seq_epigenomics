# MAST Single-Cell RNA-seq Differential Expression Analysis

This repository contains an R script for performing parallelized, block-wise differential expression analysis on single-cell RNA-seq data using the MAST package. The script includes advanced filtering steps, CPM normalization (Scanpy-compatible), gene annotation filtering, and outputs both full results and significant findings with FDR correction.

## Features

- Automated installation and loading of required CRAN and Bioconductor packages
- CPM normalization with natural log transformation for Scanpy compatibility
- Gene filtering based on single-cell expression in samples
- Removal of unwanted gene types using annotation file
- Parallelized computation over gene blocks for speed (using BiocParallel)
- Outputs: summary CSV, R objects, session info

## Requirements

- R (version >= 4.0 recommended)
- [MAST](https://bioconductor.org/packages/MAST/)
- [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/)
- [BiocParallel](https://bioconductor.org/packages/BiocParallel/)
- CRAN: `tidyverse`, `data.table`
- Input files described below

> **Note:** The script will automatically install missing packages if run interactively.

## Data Organization

Place your data in a directory structure like this *(edit to your needs)*:

your/directory/data/
    kat7_filtered_metadata.csv
    kat7_filtered_counts_matrix.csv
    MGIBatchReport_20250815_220553.txt
your/directory/output/
    (output files will be created here)


## Basic Usage 
Prepare your data files in the appropriate directories as described above.

Run the analysis:
```bash
Rscript mast_de_analysis.R
```
## Expected results:
Differential expression results: sgKat7_MAST_DE_results.csv
R objects for reproducibility: Kat7_MAST_analysis_objects_para.RData
R session info: sgKat7_session_info_para.txt

## Customization
Parallelization: The default is 36 cores. Adjust N_CORES in the script to fit your machine.

## Citation
If you use this script or adapt it, please cite:

Tang, Y.J. et al. 2025

and the MAST package:

Finak, G. et al. (2015). MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biology.

Author
Yuning J. Tang
yjtang@stanford.edu