#!/bin/bash

# Load required modules (customize or remove if not using module system)
module purge
module load fastqc

# ==== USER: Set the base data directory and subfolder containing FASTQ files ====

# Base directory for your RNA-seq project
base_dir="/your/path/to/rna_seq/normal_at2"

# Directory containing FASTQ files (customize as needed)
fastq_folder="${base_dir}/fastq/your_run_folder"

# Output directory for FastQC reports
output_dir="${base_dir}/fastqc"

# Create the output directory if it does not exist
mkdir -p "$output_dir"

# Create a list of all the FastQ files to process
fastq_files=(${fastq_folder}/*.fastq.gz)

# Check if any FastQ files were found
if [ ${#fastq_files[@]} -eq 0 ]; then
    echo "No FastQ files found in the specified directory ($fastq_folder)."
    exit 1
fi

# Process each FastQ file
for file in "${fastq_files[@]}"; do
    echo "Found file: $file"
    echo "Running FastQC on: $file"
    fastqc -o "$output_dir" "$file" & # Run FastQC in the background
done

# Wait for all background processes to finish
wait

echo "FastQC analysis completed. Reports are located in: $output_dir"