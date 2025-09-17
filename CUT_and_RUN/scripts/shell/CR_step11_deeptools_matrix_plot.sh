#!/bin/bash

###############################################################################
# Script: CR_step11_deeptools_matrix_plot.sh
# Purpose: Run deepTools computeMatrix, plotHeatmap, and plotProfile on a batch
#          of bigWig files for visualization of genome-wide signal around regions.
#
# Usage:
#   1. Edit paths, bigWig source, region annotation file as needed.
#   2. Run: bash CR_step11_deeptools_matrix_plot.sh
###############################################################################

# Load deepTools if required by HPC environment
# module purge
# module load deeptools

# ----------- CONFIGURABLE PATHS AND PARAMETERS -----------
projPath="/path/to/your_project"
bigwigDir="$projPath/bigwig"
matrixDir="$projPath/data/matrix"
plotDir="$projPath/data/plots"
cores=24    # Number of threads for computeMatrix

regionFile="/path/to/regions.tsv"           # BED/TSV of regions (TSS, genes, etc)
blacklistFile="/path/to/genome_blacklist.bed"   # Blacklist BED (optional)

mkdir -p "$matrixDir" "$plotDir"

# ----------- INPUT CHECKING -----------
if [[ ! -d "$bigwigDir" ]]; then
    echo "Error: Directory $bigwigDir does not exist."
    exit 1
fi

bigwig_files=("$bigwigDir"/*.bw)
if [[ ${#bigwig_files[@]} -eq 0 ]]; then
    echo "Error: No bigWig files found in $bigwigDir."
    exit 1
fi

for file in "${bigwig_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "Error: bigWig file not found: $file"
        exit 1
    fi
done

echo "Found bigWig files:"
for file in "${bigwig_files[@]}"; do
    echo "  $(basename "$file")"
done

if [[ ! -f "$regionFile" ]]; then
    echo "Error: Region annotation file $regionFile not found."
    exit 1
fi

# Blacklist optional: Only add if file exists
blacklist_arg=()
if [[ -f "$blacklistFile" ]]; then
    blacklist_arg=(--bl "$blacklistFile")
fi

# ----------- RUN COMPUTEMATRIX AND PLOTS -----------
echo "Running computeMatrix scale-regions ..."
computeMatrix scale-regions \
    -S "${bigwig_files[@]}" \
    -R "$regionFile" \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 0 \
    --afterRegionStartLength 3000 \
    --skipZeros \
    "${blacklist_arg[@]}" \
    -o "$matrixDir/matrix_target.mat.gz" \
    -p "$cores"

echo "Plotting heatmap ..."
plotHeatmap \
    -m "$matrixDir/matrix_target.mat.gz" \
    -out "$plotDir/heatmap.svg" \
    --sortUsing sum \
    --colorMap Reds Oranges Greys Blues Greens

echo "Plotting profile ..."
plotProfile \
    -m "$matrixDir/matrix_target.mat.gz" \
    -out "$plotDir/profile.svg" \
    --regionsLabel "" \
    --perGroup

echo "DeepTools matrix and plots complete."