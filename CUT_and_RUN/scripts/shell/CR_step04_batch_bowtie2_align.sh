#!/bin/bash

###############################################################################
# Script: CR_step04_batch_bowtie2_align.sh
# Purpose: Batch-align trimmed paired-end FASTQ files with Bowtie2
# Usage: 
#   1. Edit ref, projPath, fqPath, core count, and sample range below.
#   2. Run with: bash CR_step04_batch_bowtie2_align.sh
###############################################################################

# Optional: Load Bowtie2 module for HPC/cluster
# module purge
# module load bowtie2/2.2.0

set -e  # Exit immediately if a command fails

cores=24
ref="/path/to/Bowtie2_index/mm10"       # Bowtie2 index prefix (not .bt2 files)
projPath="/path/to/project"             # Project root path
fqPath="/path/to/trimmed_reads"         # Directory with trimmed paired reads

start_sample=1   # Index of first sample
end_sample=8     # Index of last sample

# Build sample name array
samples=()
for i in $(seq $start_sample $end_sample); do
    samples+=("CR_$i")
done

# Create output directories if missing
mkdir -p "${projPath}/alignment/sam"
mkdir -p "${projPath}/alignment/sam/bowtie2_summary"

# Alignment loop
for sample in "${samples[@]}"; do
    echo "$(date): Aligning $sample..."

    fq1="${fqPath}/merged_${sample}_1.paired.trim.fq.gz"
    fq2="${fqPath}/merged_${sample}_2.paired.trim.fq.gz"
    samOut="${projPath}/alignment/sam/${sample}_bowtie2.sam"
    summaryOut="${projPath}/alignment/sam/bowtie2_summary/${sample}_bowtie2.txt"

    # Check for input files
    if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
        echo "Error: Missing input for $sample: $fq1 or $fq2" >&2
        continue
    fi

    # Run Bowtie2
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -p "$cores" -x "$ref" -1 "$fq1" -2 "$fq2" -S "$samOut" \
        &> "$summaryOut"

    # Check exit code
    if [[ $? -ne 0 ]]; then
        echo "Error occurred during Bowtie2 alignment for sample $sample" >&2
    else
        echo "Bowtie2 alignment for sample $sample completed successfully."
    fi
done

echo "All alignments complete."