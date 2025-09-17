#!/bin/bash

###############################################################################
# Script: CR_step12_consensus_peaks.sh
# Purpose: For each target (e.g., H3K14ac, MLL1, KAT7), produce a consensus
#          BED file by intersecting replicate peak BEDs (SEACR output).
#
# Usage:
#   1. Edit projPath, targets, and replicate file names as needed.
#   2. Run: bash CR_step12_consensus_peaks.sh
###############################################################################

# module purge
# module load bedtools

projPath="/path/to/peakCalling/SEACR"
outputPath="${projPath}/intersect"
mkdir -p "$outputPath"

# -------- Example target definitions: Edit as needed -------- #
declare -A target_replicates

target_replicates=(
    [H3K14ac]="${projPath}/H3K14ac_rep1_seacr_control.peaks.stringent.bed,${projPath}/H3K14ac_rep2_seacr_control.peaks.stringent.bed"
    [MLL1]="${projPath}/MLL1_rep1_seacr_control.peaks.stringent.bed,${projPath}/MLL1_rep2_seacr_control.peaks.stringent.bed"
    [KAT7]="${projPath}/KAT7_rep1_seacr_control.peaks.stringent.bed,${projPath}/KAT7_rep2_seacr_control.peaks.stringent.bed"
)

for target in "${!target_replicates[@]}"; do
    # Parse replicate file paths
    IFS=',' read -r rep1 rep2 <<< "${target_replicates[$target]}"
    consensus_file="${outputPath}/consensus_${target}.bed"

    # Check both input files exist
    if [[ ! -f "$rep1" ]]; then
        echo "Error: Missing BED file $rep1 for $target."
        continue
    fi
    if [[ ! -f "$rep2" ]]; then
        echo "Error: Missing BED file $rep2 for $target."
        continue
    fi

    # Run bedtools intersect
    bedtools intersect -a "$rep1" -b "$rep2" > "$consensus_file"
    echo "Consensus peaks for $target written to $consensus_file"
done

echo "All consensus peak intersections completed."