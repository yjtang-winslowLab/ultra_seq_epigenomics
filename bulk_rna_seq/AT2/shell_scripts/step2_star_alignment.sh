#!/bin/bash

# Load STAR module (customize/remove as needed)
module purge
module load star

# ==== USER: Update these directories for your environment ====

# Base directory for your RNA-seq project
base_dir="/your/path/to/rna_seq/normal_at2"

# Directory containing paired FASTQ files
fastq_dir="${base_dir}/fastq/your_run_folder"

# Output directory for STAR alignment results
output_dir="${base_dir}/alignment/star_align"

# Path to STAR genome index
index="/your/path/to/Reference_genomes/star_index_mm10"

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop through all R1 FASTQ files and align paired-end reads
for fq1 in "${fastq_dir}"/*_R1_001.fastq.gz; do
    # Extract sample name (removes directory path and suffix)
    sample=$(basename "${fq1}" "_R1_001.fastq.gz")

    # Define paired-end FASTQ file (R2)
    fq2="${fastq_dir}/${sample}_R2_001.fastq.gz"

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