#!/bin/bash

###############################################################################
# Script: CR_step14_make_scaled_smoothed_bigwig.sh
# Purpose: For each sample, generate a scaled and smoothed bigWig using deepTools
# Usage: Edit paths, samples, scale factors. 
# Run: bash CR_step14_make_scaled_smoothed_bigwig.sh
###############################################################################

# module purge
# module load deeptools
# module load samtools

projPath="/path/to/project"
bamPath="${projPath}/bam"
outputPath="${projPath}/bigwigs_scaled_smoothed"
logDir="${projPath}/logs/bamCoverage_scaled_smoothed"
smooth=100           # Smoothing window size (bp)
threads=16           # Number of CPU threads to use

mkdir -p "$outputPath"
mkdir -p "$logDir"

# ========== CONFIGURE SAMPLES, TARGETS, and SCALE FACTORS ==========
samples=(Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8)
targets=(IgG H3K4me3 H3K14ac_Rep1 H3K14ac_Rep2 MLL1_Rep1 MLL1_Rep2 KAT7_Rep1 KAT7_Rep2)
scale_factors=(0.2 0.5 1.4 0.9 0.07 0.06 0.07 0.08)     # Replace with your scaling values (R/deseq/edgeR output)

# ========== MAIN PROCESSING LOOP ==========
for i in "${!samples[@]}"; do
    sample="${samples[$i]}"
    target="${targets[$i]}"
    scale_factor="${scale_factors[$i]}"

    input_bam="${bamPath}/${sample}.sorted.bam"
    output_bw="${outputPath}/${sample}_${target}_smoothed${smooth}bp_scaled.bw"
    logFile="${logDir}/${sample}_${target}_bamCoverage.log"

    if [[ ! -f "$input_bam" ]]; then
        echo "Warning: BAM not found for $sample ($target): $input_bam"
        continue
    fi

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Processing $sample ($target)..."

    bamCoverage -b "$input_bam" \
        -o "$output_bw" \
        --scaleFactor "$scale_factor" \
        --binSize 10 \
        --smoothLength "$smooth" \
        --numberOfProcessors "$threads" \
        --normalizeUsing None \
        --ignoreDuplicates \
        --extendReads \
        --centerReads &> "$logFile"

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Finished $sample ($target): $output_bw"
done

echo "All samples complete. Smoothed, scaled bigWigs are at: $outputPath"