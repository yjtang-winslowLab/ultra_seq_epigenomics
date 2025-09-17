#!/bin/bash

###############################################################################
# Script: CR_step13_annotate_consensus_peaks_homer.sh
# Purpose: Annotate consensus peak BED files using HOMER's annotatePeaks.pl.
#
# Usage:
#   - Edit projPath, input directory, genome build as needed.
#   - Run: bash CR_step13_annotate_consensus_peaks_homer.sh
###############################################################################

# module purge
# module load homer/4.11

# Confirm HOMER annotatePeaks.pl is available
if ! command -v annotatePeaks.pl &> /dev/null; then
    echo "Error: HOMER module not loaded or annotatePeaks.pl not in PATH. Exiting."
    exit 1
fi

projPath="/path/to/project/peakCalling/SEACR"
inputBedDir="${projPath}/processed_bedfiles"
outputPath="${projPath}/annotation"
mkdir -p "$outputPath"

# Collect BED files to annotate
inputFiles=()
for f in "$inputBedDir"/*.bed; do
    if [[ -f "$f" ]]; then
        inputFiles+=("$f")
    fi
done

if [[ ${#inputFiles[@]} -eq 0 ]]; then
    echo "No BED files found in $inputBedDir."
    exit 1
fi

echo "BED files to annotate:"
for f in "${inputFiles[@]}"; do
    echo "  $f"
done

# Genome build
genome="mm10"

for bedFile in "${inputFiles[@]}"; do
    base="$(basename "$bedFile" .bed)"
    outputFile="$outputPath/${base}_annotated.txt"
    annStatsFile="$outputPath/${base}_annotation_stats.txt"

    echo "Annotating $bedFile ..."
    annotatePeaks.pl "$bedFile" "$genome" -annStats "$annStatsFile" -p > "$outputFile"
    if [ $? -eq 0 ]; then
        echo "  Annotated peaks written: $outputFile"
        echo "  Annotation stats written: $annStatsFile"
    else
        echo "Error: Failed to annotate $bedFile. Check input and HOMER installation."
    fi
done

echo "All annotation jobs complete."