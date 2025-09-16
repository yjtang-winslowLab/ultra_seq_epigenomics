#!/bin/bash

# Script to format merged BED files into SAF and BED formats with peak IDs
# for downstream analysis (e.g., featureCounts, differential analysis).
# Requires bedops module loaded.

# Load required modules (adjust according to your HPC environment)
module purge
module load bedops/2.4.41

# Define base directory (edit this path as needed)
BASE_DIR="path/to/your/output_directory"

# Define input and output files
INPUT_BED="${BASE_DIR}/All_Samples.fwp.filter.non_overlapping.bed"
OUTPUT_SAF="${BASE_DIR}/All_Samples.fwp.filter.non_overlapping.saf"
OUTPUT_PEAKID_BED="${BASE_DIR}/All_Samples.fwp.filter.non_overlapping.peakIds.bed"

# Convert BED file to SAF format (for featureCounts)
awk -F $'\t' 'BEGIN {OFS = FS; nr=0} { 
    $2=$2+1; 
    peakid="all_samples_summits_"++nr;  
    print peakid,$1,$2,$3,"."
}' "$INPUT_BED" > "$OUTPUT_SAF"

# Convert BED file to BED format with peak IDs (for downstream analysis)
awk -F $'\t' 'BEGIN {OFS = FS; nr=0} { 
    $2=$2+1; 
    peakid="all_samples_summits_"++nr;  
    print $1,$2,$3,peakid,"0","."
}' "$INPUT_BED" > "$OUTPUT_PEAKID_BED"
