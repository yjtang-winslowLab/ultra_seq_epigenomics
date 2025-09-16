#!/bin/bash

# Load required modules
module purge
module load fastqc

# ==== USER: Set these paths to your own data location ====

# Define the base directory for input files
base_dir="/your/reference/data/rna_seq"

# Define the output directory for FastQC reports
output_dir="${base_dir}/fastqc"

# Create the output directory if it does not exist
mkdir -p "$output_dir"

# Create a list of all the FastQ files to process (use .fq.gz extension)
fastq_files=(${base_dir}/*.fq.gz)

# Check if any FastQ files were found
if [ ${#fastq_files[@]} -eq 0 ]; then
    echo "No FastQ files found in the specified directory: $base_dir"
    exit 1
fi

# Process each FastQ file
for file in "${fastq_files[@]}"; do
    echo "Found file: $file"
    echo "Running FastQC on: $file"
    fastqc -o "$output_dir" "$file" &
done

# Wait for all background processes to finish
wait

echo "FastQC analysis completed. Reports are located in: $output_dir"