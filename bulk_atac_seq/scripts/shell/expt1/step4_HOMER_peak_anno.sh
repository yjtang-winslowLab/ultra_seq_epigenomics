#!/bin/sh

# Script to annotate peaks using HOMER's annotatePeaks.pl
# Requires HOMER module loaded (annotatePeaks.pl)
# Adjust paths, filenames, and genome according to your environment

# Load required modules (adjust according to your HPC environment)
module purge
module load homer/4.11

# Define genome (edit as needed, e.g., mm10, hg38)
GENOME="mm10"

# Define input and output directories (edit these paths as needed)
INPUT_DIR="path/to/your/input_directory/"
OUTPUT_DIR="path/to/your/output_directory/"

# List of datasets to annotate (edit this array with your filenames, without extensions)
datasets=("peaks_for_HOMER") # single file here, adjust as needed

# Loop through each dataset and run HOMER annotation
for dataset in "${datasets[@]}"
do
    echo "Running HOMER annotatePeaks on ${dataset}"

    # Define input and output filenames
    INPUT_FILE="${INPUT_DIR}${dataset}.txt"
    OUTPUT_FILE="${OUTPUT_DIR}${dataset}_annotated_peaks.txt"

    # Run HOMER annotatePeaks.pl
    annotatePeaks.pl "$INPUT_FILE" "$GENOME" > "$OUTPUT_FILE"

    echo "Annotation complete for ${dataset}. Results saved to ${OUTPUT_FILE}"
done