#!/bin/bash

###############################################################################
# Script: CR_step06_sam_to_fragments.sh
# Purpose: For each sample: convert mapped SAM to BAM, then BEDPE, then to sorted
#          single BED format for fragment-level coverage/scanning.
#
# Usage: 
#   1. Edit projPath below. Make sure SAM files are present.
#   2. Run: bash CR_step06_sam_to_fragments.sh
###############################################################################

# Optional: load modules for cluster
# module purge
# module load samtools
# module load bedtools

set -e  # Stop on first error

projPath="/path/to/project"
bamPath="${projPath}/bam"
bedPath="${projPath}/bed"
samPath="${projPath}/alignment/sam"

# Create output directories if needed
mkdir -p "$bamPath" "$bedPath" "$samPath"

start_sample=1
end_sample=8      # Adjust for your sample count
samples=()
for i in $(seq $start_sample $end_sample); do
    samples+=("CR6_$i")
done

for sample in "${samples[@]}"; do
    echo "Processing $sample..."

    samFile="${samPath}/${sample}_bowtie2.sam"
    bamFile="${bamPath}/${sample}_bowtie2.mapped.bam"
    bedpeFile="${bedPath}/${sample}_bowtie2.bed"
    fragFile="${bedPath}/${sample}_bowtie2.fragments.bed"

    # Check for SAM file
    if [[ ! -f "$samFile" ]]; then
        echo "Error: SAM file for $sample not found: $samFile" >&2
        continue
    fi

    # 1. Filter mapped read pairs to BAM
    samtools view -bS -F 0x04 "$samFile" > "$bamFile"

    # 2. Convert BAM to BEDPE (paired-end)
    bedtools bamtobed -i "$bamFile" -bedpe > "$bedpeFile"

    # 3. Keep read pairs on the same chromosome and <1000bp, extract start/end
    awk -v OFS="\t" '$1==$4 && $6-$2 < 1000 {print $1, $2, $6}' "$bedpeFile" | \
        sort -k1,1 -k2,2n > "$fragFile"

    echo "Generated fragment BED: $fragFile"
done

echo "All SAM → BAM → BED fragment processing complete."