#!/bin/bash

###############################################################################
# Script: CR_step10_bam_to_bigwig.sh
# Purpose: Sort and index BAM files, then generate bigWig files with deepTools.
#
# Usage:
#   1. Edit project paths and sample_map for your dataset.
#   2. Run: bash CR_step10_bam_to_bigwig.sh
###############################################################################

# module purge
# module load samtools
# module load deeptools

projPath="/path/to/your_project"
bamDir="${projPath}/bam"
outputDir="${projPath}/bigwig"
mkdir -p "$outputDir"

# Sample mapping: sample ID → display name in output file
declare -A sample_map
sample_map=(
    [Sample1]="IgG"
    [Sample2]="H3K4me3"
    [Sample3]="H3K14ac_Rep1"
    [Sample4]="H3K14ac_Rep2"
    [Sample5]="MLL1_Rep1"
    [Sample6]="MLL1_Rep2"
    [Sample7]="KAT7_Rep1"
    [Sample8]="KAT7_Rep2"
)

for sample in "${!sample_map[@]}"; do
    histName="${sample_map[$sample]}"
    input_bam="${bamDir}/${sample}_mapped.bam"
    sorted_bam="${bamDir}/${sample}.sorted.bam"
    output_bw="${outputDir}/${histName}_raw.bw"

    if [[ ! -f "$input_bam" ]]; then
        echo "Warning: BAM file not found for $sample: $input_bam"
        continue
    fi

    echo "Sorting $input_bam..."
    samtools sort -o "$sorted_bam" "$input_bam"

    echo "Indexing $sorted_bam..."
    samtools index "$sorted_bam"

    echo "Generating bigWig for $histName..."
    bamCoverage -b "$sorted_bam" -o "$output_bw"

    echo "Done: $output_bw"
done

echo "All BAMs converted to bigWig files in $outputDir."