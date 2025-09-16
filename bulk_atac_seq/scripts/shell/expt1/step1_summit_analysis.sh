#!/bin/sh

# Script to perform iterative overlap peak merging for ATAC-seq data
# Uses createIterativeOverlapPeakSet.R from Corces Lab:
# https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging

# Load required modules (adjust according to your HPC environment)
module purge
module load R/4.2.2

# Define variables (edit these paths and parameters as needed)
METADATA="path/to/your/metadata.txt"                 # Metadata file describing samples
MACS2_DIR="path/to/your/macs2_peaks_directory"       # Directory containing MACS2 peak files
OUTPUT_DIR="path/to/your/output_directory"           # Directory to store output files
SUFFIX="_summits.bed"                                # Suffix of MACS2 summit files
BLACKLIST="path/to/your/genome_blacklist.bed"        # Genome blacklist regions to exclude
GENOME="mm10"                                        # Reference genome (e.g., mm10, hg38)
SUMMITS_PER_MERGE=5                                  # Number of summits per merge iteration
MERGE_RULE="2"                                       # Merge rule (e.g., minimum number of overlaps)
EXTEND_SIZE=250                                      # Extend summits by this many base pairs

# Run the iterative overlap peak merging R script
Rscript scripts/createIterativeOverlapPeakSet.R \
  --metadata "$METADATA" \
  --macs2dir "$MACS2_DIR" \
  --outdir "$OUTPUT_DIR" \
  --suffix "$SUFFIX" \
  --blacklist "$BLACKLIST" \
  --genome "$GENOME" \
  --spm "$SUMMITS_PER_MERGE" \
  --rule "$MERGE_RULE" \
  --extend "$EXTEND_SIZE"
