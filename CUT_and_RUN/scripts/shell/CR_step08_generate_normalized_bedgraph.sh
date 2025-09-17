#!/bin/bash

###############################################################################
# Script: CR_step08_generate_normalized_bedgraph.sh
# Purpose: Normalize CUT&RUN (or ChIP-seq, ATAC-seq) .bed fragment files using 
#          sample-specific scale factors and generate bedGraph coverage files.
# Requirements:
#   - bedtools (for genomecov)
#   - bash >= 4
# Usage:
#   Edit the paths and sample info in the configuration section below.
###############################################################################

# Load modules if required (for cluster environments, otherwise comment/remove)
# module purge
# module load bedtools

###########################
# Configuration - EDIT ME #
###########################

projPath="/path/to/project"
inputPath="$projPath/analysis/bed"
chromSize="/path/to/reference.genome.chrom.sizes"

outputDir="$projPath/analysis/bedgraph"
mkdir -p "$outputDir"

# Scale factors (replace with your computed values)
# This is calculated in R. Check the other folder for details. 
scale_factors=(0.10 0.25 1.00 0.55 0.07 0.60 0.88 0.95)

# Sample identifiers (edit to match your sample names/IDs)
samples=("SampleA" "SampleB" "SampleC" "SampleD" "SampleE" "SampleF" "SampleG" "SampleH")

# Target names (edit to match what you want in the output filenames)
targets=("Control" "H3K4me3" "H3K14ac_1" "H3K14ac_2" "MLL1_1" "MLL1_2" "KAT7_1" "KAT7_2")

###################################
# Processing - Do not edit below  #
###################################

for index in "${!samples[@]}"; do
    sample="${samples[$index]}"
    target="${targets[$index]}"
    scale_factor="${scale_factors[$index]}"
    
    # Input fragment .bed file for this sample
    inputFile="$inputPath/${sample}_fragments.bed"
    outputFile="$outputDir/${target}_fragments.normalized.bedgraph"
    
    # Check input file exists
    if [[ ! -f "$inputFile" ]]; then
        echo "Warning: Input file $inputFile not found, skipping ${sample}"
        continue
    fi
    
    # Generate normalized bedGraph
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Generating normalized bedGraph for ${sample} (${target}), scale=${scale_factor}"
    bedtools genomecov -bg -scale "$scale_factor" -i "$inputFile" -g "$chromSize" > "$outputFile"
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Finished ${sample} (${target})"
done

echo "All samples processed. Output in $outputDir"