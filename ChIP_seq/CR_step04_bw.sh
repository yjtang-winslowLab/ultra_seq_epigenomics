#!/bin/bash

# Load necessary modules
module purge
module load deeptools
module load samtools

# Define paths -- Make sure to change these!
projPath="/your/directory/Keuh"
bamPath="${projPath}/keuh_alignment/bam"
bwPath="${projPath}/keuh_bw"
logPath="${projPath}/logs/bam_to_bw_logs"

# Create output directories
mkdir -p "$bwPath"
mkdir -p "$logPath"

# Effective genome size for mouse mm10
effectiveGenomeSize=2652783500

# Fragment length (average fragment size). Adjust if known from your data.
fragmentLength=100

# Define samples and labels
declare -A samples=(
    ["SRR22028893"]="igG"
    ["SRR22028897"]="H3K14ac_rep4"
    ["SRR22028898"]="H3K14ac_rep3"
    ["SRR22028899"]="H3K14ac_rep1"
    ["SRR22028900"]="H3K14ac_rep2"
)

for sample in "${!samples[@]}"; do
    sorted_bam="${bamPath}/${sample}.sorted.bam"
    bwFile="${bwPath}/${samples[$sample]}_RPGC_norm.bw"
    logFile="${logPath}/${samples[$sample]}_bamCoverage.log"

    {
        echo "==============================================="
        echo "$(date): Starting BigWig generation for ${samples[$sample]}"
        echo "Sample ID: ${sample}"
        echo "Sorted BAM: ${sorted_bam}"
        echo "Output BigWig: ${bwFile}"
        echo "==============================================="

        # Check if BAM file exists
        if [[ ! -f "$sorted_bam" ]]; then
            echo "Error: Sorted BAM file $sorted_bam not found!"
            continue
        fi

        # Index BAM if needed
        if [[ ! -f "${sorted_bam}.bai" ]]; then
            echo "$(date): Indexing BAM file..."
            samtools index "$sorted_bam"
        else
            echo "$(date): BAM index already exists."
        fi

        echo "$(date): Running bamCoverage..."
        bamCoverage \
            --bam "$sorted_bam" \
            --outFileName "$bwFile" \
            --binSize 10 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize "$effectiveGenomeSize" \
            --extendReads "$fragmentLength" \
            --numberOfProcessors 8

        if [[ $? -eq 0 ]]; then
            echo "$(date): Successfully generated BigWig for ${samples[$sample]}."
        else
            echo "$(date): Error occurred during bamCoverage for ${samples[$sample]}!"
        fi

        echo "==============================================="
        echo ""
    } &> "$logFile"

done

echo "$(date): All normalized BigWig files generated. Logs are stored in ${logPath}."