#!/bin/bash

# Load required module (customize/remove as needed for your environment)
module purge
module load subread

# ==== USER: Set these to your project directory and annotation ====

# Base directory for your RNA-seq project
base_dir="/your/path/to/rna_seq/normal_at2"

# Directory containing BAM files from STAR alignment
alignment_dir="${base_dir}/alignment/star_align"

# Output directory for featureCounts results
output_dir="${base_dir}/counts"

# Path to annotation GTF file
gtf_file="/your/path/to/Reference_genomes/gencode.vM10.annotation.gtf"

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

# Run featureCounts on all BAM files in the alignment directory
featureCounts -T 16 -p \
  -t exon \
  -g gene_id \
  -a "${gtf_file}" \
  -o "${output_dir}/featureCounts_results.txt" \
  "${alignment_dir}"/*.bam