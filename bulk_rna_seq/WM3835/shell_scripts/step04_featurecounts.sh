#!/bin/bash

# Load required modules
module purge
module load subread

# ==== USER: Set these directories and filenames for your data/reference ====

# Base directory for input files
base_dir="/your/reference/data/rna_seq"

# Directory containing alignment BAM files
alignment_dir="${base_dir}/alignment/star_align"

# Output directory for featureCounts results
output_dir="${base_dir}/counts"

# Path to annotation GTF file
gtf_file="/your/reference/genomes/gencode.vM10.annotation.gtf"

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

# List of BAM files to count - update with your own filenames
bam_files=(
  "KP787_3_D1_Aligned.sortedByCoord.out.bam"
  "KP787_3_D3_Aligned.sortedByCoord.out.bam"
  "KP787_3_WM1_Aligned.sortedByCoord.out.bam"
  "KP787_3_WM3_Aligned.sortedByCoord.out.bam"
)

# Prepend alignment directory path to each BAM file
bam_paths=()
for bam in "${bam_files[@]}"; do
  bam_paths+=("${alignment_dir}/${bam}")
done

# Run featureCounts in parallel mode (-T), for paired-end data (-p)
featureCounts -T 16 -p \
  -t exon \
  -g gene_id \
  -a "${gtf_file}" \
  -o "${output_dir}/featureCounts_results.txt" \
  "${bam_paths[@]}"