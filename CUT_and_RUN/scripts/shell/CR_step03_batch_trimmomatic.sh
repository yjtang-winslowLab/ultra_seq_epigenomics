#!/bin/bash

###############################################################################
# Script: CR_step03_batch_trimmomatic.sh
# Purpose: Batch trimming of paired-end FASTQ files with Trimmomatic
# Usage:
#   1. Set projPath (input), outputDir (output), and adapterPath below.
#   2. Adjust sample range (start_sample, end_sample), threads, quality, and minLength as needed.
#   3. Run: bash CR_step03_batch_trimmomatic.sh
###############################################################################

# Optional: load Trimmomatic module for cluster environments
# module purge
# module load trimmomatic/0.39

# ----- CONFIGURE THESE PATHS -----
projPath="/path/to/raw_data"              # Parent directory with sample folders
outputDir="/path/to/trimmed_data"         # Output directory for trimmed FASTQ files
adapterPath="/path/to/TruSeq2-PE.fa"      # Path to adapters file (adjust as needed)
trimmomaticJar="/path/to/trimmomatic-0.39.jar"  # Path to Trimmomatic jar

threads=24        # Number of threads
quality=20        # Sliding window quality cutoff
minLength=50      # Minimum read length after trimming
start_sample=1    # Index of first sample
end_sample=8      # Index of last sample

mkdir -p "$outputDir"

# Array of sample names
samples=()
for i in $(seq $start_sample $end_sample); do
    samples+=("CR_$i")
done

# -------- Function: run Trimmomatic on a sample ----------
run_trimmomatic() {
    local sample=$1
    echo "$(date): Processing $sample..."

    # Input FASTQ files
    in1="${projPath}/${sample}/merged_${sample}_1.fq.gz"
    in2="${projPath}/${sample}/merged_${sample}_2.fq.gz"

    # Check input files exist
    if [[ ! -f "$in1" || ! -f "$in2" ]]; then
        echo "Error: Missing input for $sample: $in1 or $in2"
        return 1
    fi

    # Output files
    out1_paired="${outputDir}/merged_${sample}_1.paired.trim.fq.gz"
    out1_unpaired="${outputDir}/merged_${sample}_1.unpaired.trim.fq.gz"
    out2_paired="${outputDir}/merged_${sample}_2.paired.trim.fq.gz"
    out2_unpaired="${outputDir}/merged_${sample}_2.unpaired.trim.fq.gz"

    # Run Trimmomatic
    java -jar "${trimmomaticJar}" PE \
        -threads "$threads" -phred33 \
        "$in1" "$in2" \
        "$out1_paired" "$out1_unpaired" \
        "$out2_paired" "$out2_unpaired" \
        ILLUMINACLIP:"$adapterPath":2:40:15 \
        SLIDINGWINDOW:4:"$quality" \
        MINLEN:"$minLength"

    # Check for success
    if [[ $? -ne 0 ]]; then
        echo "Error: Trimmomatic failed for $sample!"
        return 1
    fi

    echo "Trimming complete for ${sample} at $(date)."
}

# -------- Main Loop --------
for sample in "${samples[@]}"; do
    run_trimmomatic "$sample"
done

echo "All processing complete."