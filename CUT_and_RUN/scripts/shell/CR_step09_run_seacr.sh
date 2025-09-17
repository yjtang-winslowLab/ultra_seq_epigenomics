#!/bin/bash

###############################################################################
# Script: run_seacr.sh
# Purpose: Run SEACR for peak calling on normalized bedgraph files for multiple targets, 
#          using a control bedgraph file.
#
# Usage:
#   1. Edit projPath, seacrScript, targets, and control sample as needed.
#   2. Run: bash run_seacr.sh
#
###############################################################################

# Load necessary modules (adapt/remove if not using modules or cluster)
# module purge
# module load SEACR/1.3
# module load bedtools
# module load R

# ---- CONFIGURE THESE VARIABLES FOR YOUR PROJECT ----
projPath="/path/to/project"
seacrScript="/path/to/SEACR_1.3.sh"

# List of target samples (edit as needed)
targets=(
    "Sample1"
    "Sample2"
    "Sample3"
    "Sample4"
    "Sample5"
    "Sample6"
    "Sample7"
)

control="ControlSample"   # Control name in bedgraph file (e.g., IgG, Input)

inputDir="$projPath/bedgraph"
outputDir="$projPath/peakCalling/SEACR"
mkdir -p "$outputDir"

for target in "${targets[@]}"; do
    inputBedgraph="${inputDir}/${target}_fragments.normalized.bedgraph"
    controlBedgraph="${inputDir}/${control}_fragments.normalized.bedgraph"

    # Check files exist before running
    if [[ ! -f "$inputBedgraph" ]]; then
        echo "Warning: $inputBedgraph does not exist. Skipping $target."
        continue
    fi
    if [[ ! -f "$controlBedgraph" ]]; then
        echo "Warning: $controlBedgraph does not exist. Skipping $target."
        continue
    fi

    # Run SEACR for peaks vs control
    bash "$seacrScript" "$inputBedgraph" "$controlBedgraph" non stringent \
        "$outputDir/${target}_seacr_control.peaks"

    # Optionally: Run SEACR for top N% threshold peaks (uncomment if needed)
    # bash "$seacrScript" "$inputBedgraph" 0.01 non stringent \
    #    "$outputDir/${target}_seacr_top0.01.peaks"
    
    echo "SEACR peak calling complete for $target."
done

echo "All SEACR peak calling jobs finished."