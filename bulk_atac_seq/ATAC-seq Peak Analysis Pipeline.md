ATAC-seq Peak Analysis Pipeline

Author: Y.J. Tang, Winslow Lab, Stanford University

This repository contains scripts for processing and analyzing bulk ATAC-seq data described in Tang et al. (Figure 5 and associated Supplemental Figures). The ATAC-seq experiments consist of two separate experiments:
Expt1: Includes sgRNAs targeting MLL1 complex members and Meaf6.
Expt2: Includes sgRNAs targeting Kat7.

Preprocessing of FASTQ files and peak-calling were performed by Samuel H. Kim (MD/Phd., Will Greenleaf Lab, Stanford University) based on scripts from Georgi Marinov's GitHub repository.

Please refer to the provided metadata file for detailed sample information.

Overview
The pipeline includes the following shell scripts:
step1_summit_analysis.sh: Performs iterative overlap peak merging using an R script from the Corces Lab.
step2_format_bedfiles.sh: Formats merged BED files into SAF and BED formats with peak IDs for downstream analysis.
step3_summits_counts.sh: Quantifies reads overlapping peaks using featureCounts.
step4_HOMER_peak_anno.sh: Annotates peaks using HOMER's annotatePeaks.pl.
step5_HOMER_MotifAnalysis.sh: Identifies transcription factor motifs enriched in peak regions. (This script is specifically included for Expt1.)
To investigate the overlap between the MLL1 complex and the HBO1 complex, the combined analysis uses peak coordinates from Expt1 to quantify peaks in Expt2.

Requirements
Ensure the following software and modules are installed and available in your HPC environment:
R (tested with version 4.2.2)
bedops (tested with version 2.4.41)
subread (featureCounts) (tested with version 2.0.3)
HOMER (tested with version 4.11)

Customization
Before running the scripts, edit the paths and parameters at the top of each script to match your computing environment and data locations.

Citation and Acknowledgments
If you use this pipeline, please cite the original manuscript (Tang et al.) and the following tools and scripts:
Iterative Overlap Peak Merging: Corces Lab GitHub
featureCounts: Subread
HOMER: HOMER website
Preprocessing and Peak-calling scripts: Georgi Marinov's GitHub repository
