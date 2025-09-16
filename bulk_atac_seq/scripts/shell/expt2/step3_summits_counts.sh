#!/bin/sh

# Load required modules (adjust according to your HPC environment)
module purge
module load subread/2.0.3

# Define directories (edit these paths as needed)
bam_base="path/to/your/bam_directory"
saf_base="path/to/your/ouput_directory"

# Run featureCounts to quantify reads overlapping peaks
featureCounts -p -F SAF \
  -a ${saf_base}All_Samples.fwp.filter.non_overlapping.saf \
  --fracOverlap 0.2 \
  -o ${saf_base}All_Samples_summit.counts \
  ${bam_base}Bulk_ATAC_21.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_22.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_23.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_24.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_25.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_26.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_27.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_28.noMT.filtered.deduped.bam
