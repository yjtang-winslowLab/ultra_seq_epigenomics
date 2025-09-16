#!/bin/sh

# Load required modules (adjust according to your HPC environment)
module purge
module load R/4.2.2

# Define variables (edit these paths and parameters as needed)
METADATA="path/to/your/metadata.txt"
MACS2_DIR="path/to/your/macs2_peaks_directory"
OUTPUT_DIR="path/to/your/output_directory"
SUFFIX="_summits.bed"
BLACKLIST="path/to/your/genome_blacklist.bed"
GENOME="mm10"
SUMMITS_PER_MERGE=5
MERGE_RULE="2"
EXTEND_SIZE=250

# Run the iterative overlap peak merging R script
# createIterativeOverlapPeakSet.R is from https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging
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