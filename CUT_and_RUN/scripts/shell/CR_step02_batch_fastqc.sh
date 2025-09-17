#!/bin/bash

###############################################################################
# Script: CR_step02_batch_fastqc.sh
# Purpose: Run FastQC in parallel on merged FASTQ files in sample directories
# Usage:
#   1. Set base_dir and output_dir below.
#   2. Run: bash CR_step02_batch_fastqc.sh
###############################################################################

# If using on cluster: Load FastQC module (remove if unnecessary)
# module purge
# module load fastqc

base_dir="/path/to/raw_data"         # Parent directory containing CR_1 ... CR_8
output_dir="/path/to/fastqc_output"  # Where FastQC reports will be stored

# Create output directory if missing
mkdir -p "$output_dir"

echo "Starting FastQC on merged FASTQ files..."

for lib in {1..8}; do
    # Loop over both _1 and _2 files for each sample
    for read_id in 1 2; do
        fq="${base_dir}/CR_${lib}/merged_CR_${lib}_${read_id}.fq.gz"
        if [[ -f "$fq" ]]; then
            echo "Found file: $fq"
            echo "Running FastQC..."
            fastqc -o "$output_dir" "$fq" &
        else
            echo "File not found: $fq"
        fi
    done
done

# Wait for all background jobs to finish
wait

echo "FastQC analysis completed. Reports are in: $output_dir"