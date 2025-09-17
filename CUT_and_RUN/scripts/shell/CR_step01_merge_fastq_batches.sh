#!/bin/bash

###############################################################################
# Script: CR_step01_merge_fastq_batches.sh
# Purpose: Merge multiple FASTQ.gz files in subfolders for batch sequencing
# Usage:
#   - Set base_dir to the parent directory containing CR_1 ... CR_8 subfolders
#   - Script merges all *_1.fq.gz into merged_CR_${i}_1.fq.gz and *_2.fq.gz
#   - Outputs status for each sample
###############################################################################

set -e  # Exit immediately if a command fails.

# Modify as appropriate:
base_dir="/path/to/raw_data_directory"

# Loop over sample directories, e.g., CR_1 ... CR_8
for i in {1..8}; do
    subdir="$base_dir/CR_$i"
    echo "Processing $subdir..."

    # Check the directory exists
    if [ ! -d "$subdir" ]; then
        echo "Directory $subdir does not exist. Skipping."
        continue
    fi

    cd "$subdir"

    # Merge forward (_1) reads
    if ls *_1.fq.gz 1> /dev/null 2>&1; then
        cat *_1.fq.gz > merged_CR_${i}_1.fq.gz
        echo "Merged *_1.fq.gz into merged_CR_${i}_1.fq.gz"
    else
        echo "No *_1.fq.gz files found in $subdir"
    fi

    # Merge reverse (_2) reads
    if ls *_2.fq.gz 1> /