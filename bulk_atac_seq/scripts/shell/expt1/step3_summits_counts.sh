#!/bin/sh

# Script to quantify reads overlapping peaks using featureCounts
# Requires subread module loaded (featureCounts)
# Adjust paths and filenames according to your environment

# Load required modules (adjust according to your HPC environment)
module purge
module load subread/2.0.3

# Define directories (edit these paths as needed)
BAM_DIR="path/to/your/bam_directory"
OUTPUT_DIR="path/to/your/output_directory"

# Define input SAF file and output counts file
SAF_FILE="${OUTPUT_DIR}/All_Samples.fwp.filter.non_overlapping.saf"
COUNTS_FILE="${OUTPUT_DIR}/All_Samples_summit.counts"

# List of BAM files to quantify (edit as needed)
BAM_FILES=(
  "${BAM_DIR}/sgMeaf6_1.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgMeaf6_2.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgMeaf6_3.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgKmt2a_1.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgKmt2a_2.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgKmt2a_3.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgPsip1_1.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgPsip1_2.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgPsip1_3.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgSafe14_1.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgSafe23_1.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgSafe23_2.noMT.filtered.deduped.bam"
  "${BAM_DIR}/sgSafe23_3.noMT.filtered.deduped.bam"
)

# Run featureCounts to quantify reads overlapping peaks
featureCounts -p -F SAF \
  -a "$SAF_FILE" \
  --fracOverlap 0.2 \
  -o "$COUNTS_FILE" \
  "${BAM_FILES[@]}"