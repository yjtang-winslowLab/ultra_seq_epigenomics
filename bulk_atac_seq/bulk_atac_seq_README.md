# ATAC-seq Peak Analysis Pipeline

**Author:** Y.J. Tang  
**Affiliation:** Winslow Lab, Dept. of Genetics, Stanford University

---

## Overview

This repository contains scripts for processing and analyzing bulk ATAC-seq data from Tang et al. (Figure 5 and associated Supplemental Figures).  
The ATAC-seq experiments consist of two separate sets:

- **Expt1:** Contains sgRNAs targeting MLL1 complex members and Meaf6.
- **Expt2:** Contains sgRNAs targeting Kat7.

Preprocessing of FASTQ files and peak-calling were performed by Samuel H. Kim (MD/PhD., Will Greenleaf Lab, Stanford University), using scripts from Georgi Marinov's GitHub repository.

See the provided metadata file for detailed sample information.

---

## Pipeline Steps

The pipeline consists of the following shell scripts:

1. **step1_summit_analysis.sh**  
   Iterative overlap peak merging using an R script from the Corces Lab.

2. **step2_format_bedfiles.sh**  
   Formats merged BED files into SAF and BED formats, assigning peak IDs for downstream analysis.

3. **step3_summits_counts.sh**  
   Quantifies reads overlapping peaks using featureCounts.

4. **step4_HOMER_peak_anno.sh**  
   Annotates peaks using HOMER's `annotatePeaks.pl`.

5. **step5_HOMER_MotifAnalysis.sh**  
   Identifies transcription factor motifs enriched in peak regions.  
   *(This script is included specifically for Expt1.)*

**Combined Analysis:**  
To investigate overlaps between the MLL1 complex and the HBO1 complex, the pipeline uses peak coordinates from Expt1 to quantify overlaps in Expt2.

---

## Requirements

Ensure the following software/tools are installed and available in your HPC environment:

- `R` (tested with v4.2.2)
- `bedops` (tested with v2.4.41)
- `subread` (`featureCounts`) (tested with v2.0.3)
- `HOMER` (tested with v4.11)

---

## Usage & Customization

Before running the scripts, **edit the paths and parameters at the top of each script** to match your computing environment and data locations.

---

## Citation & Acknowledgments

If you use this pipeline, please cite:

- Tang et al. (original manuscript)
- Tools & scripts:
  - Iterative Overlap Peak Merging: [Corces Lab GitHub](https://github.com/corceslab)
  - featureCounts: [Subread](http://subread.sourceforge.net/)
  - HOMER: [HOMER website](http://homer.ucsd.edu/homer/)
  - Preprocessing & Peak-calling scripts: [Georgi Marinov's GitHub](https://github.com/GeorgiMarinov)

---

For questions, please contact Y.J. Tang at Stanford University.