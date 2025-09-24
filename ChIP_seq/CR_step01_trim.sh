#!/bin/bash

# Requirements:
#   - Bash, Trimmomatic v0.39, module system (adjust as needed)
#   - Adjust paths as needed

set -e

# Load necessary modules
module purge
module load trimmomatic/0.39

# Define project paths and parameters
base_dir="/your/directory/fastq"
outputDir="/your/directory/trim"
adapterPath="/scg/apps/software/trimmomatic/0.39/adapters/TruSeq3-SE.fa" # single-end reads
threads=24
quality=20
minLength=50

# Create output directories
mkdir -p "$outputDir/logs"

# Define the list of FastQ files to process
files=(
    "$base_dir/SRR22028893.fastq.gz"
    "$base_dir/SRR22028897.fastq.gz"
    "$base_dir/SRR22028898.fastq.gz"
    "$base_dir/SRR22028899.fastq.gz"
    "$base_dir/SRR22028900.fastq.gz"
)

# Function to run Trimmomatic with error checking and logging
run_trimmomatic() {
    local file=$1
    local sample=$(basename "$file" .fastq.gz)

    echo "$(date): Starting trimming for $sample..."

    java -jar /scg/apps/software/trimmomatic/0.39/trimmomatic-0.39.jar SE \
        -threads "$threads" -phred33 \
        "$file" \
        "${outputDir}/${sample}.trim.fq.gz" \
        ILLUMINACLIP:"$adapterPath":2:40:15 SLIDINGWINDOW:4:"$quality" MINLEN:"$minLength" \
        > "${outputDir}/logs/${sample}_trimmomatic.log" 2>&1

    if [[ $? -ne 0 ]]; then
        echo "Error: Trimmomatic failed for sample $sample! Check log file."
        return 1
    fi

    echo "$(date): Trimming complete for $sample."
}

# Loop through files and run trimming
for file in "${files[@]}"; do
    run_trimmomatic "$file"
done

echo "$(date): All processing complete."