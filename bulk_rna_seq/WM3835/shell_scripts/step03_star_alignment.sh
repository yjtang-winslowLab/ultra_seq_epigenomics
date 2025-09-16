#!/bin/bash
# Example STAR paired-end alignment script

# Load required modules
module purge
module load star

# ==== USER: Set these directories for your own data/reference ====

# Base directory for input files
base_dir="/your/reference/data/rna_seq"

# Directory containing FASTQ files
fastq_dir="${base_dir}"

# Output directory for STAR alignment results
output_dir="${base_dir}/alignment/star_align"

# Path to the STAR genome index
index="/your/reference/genomes/star_index_mm10"

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop through all R1 FASTQ files and align paired-end reads
for fq1 in "${fastq_dir}"/*_1.fq.gz; do
    # Extract sample name by removing directory path and suffix "_1.fq.gz"
    sample=$(basename "${fq1}" "_1.fq.gz")

    # Define the paired-end FASTQ file (R2)
    fq2="${fastq_dir}/${sample}_2.fq.gz"

    # Check if paired-end file exists
    if [[ ! -f "${fq2}" ]]; then
        echo "Paired file ${fq2} not found for sample ${sample}, skipping..."
        continue
    fi

    # Echo filenames to verify correctness
    echo "Sample: ${sample}"
    echo "Read 1: ${fq1}"
    echo "Read 2: ${fq2}"
    echo "-----------------------------------"

    # Run STAR alignment
    STAR --genomeDir "${index}" \
         --runMode alignReads \
         --readFilesIn "${fq1}" "${fq2}" \
         --runThreadN 16 \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --outReadsUnmapped Fastx \
         --outFileNamePrefix "${output_dir}/${sample}_"

done