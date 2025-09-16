#!/bin/sh

# Load required modules (adjust according to your HPC environment)
module purge
module load subread/2.0.3

## Generate read counts using genomic regions (SAF file) from one batch (expt1: sgX vs sgSafe)
## applied to BAM files from another batch (expt2: sgKat7 vs sgSafe).

# Define directories (edit these paths as needed)
bam_base=/your/directory/to/bam/from/expt2
saf_base=/your/directory/to/saf/from/expt1 # SAF file from expt1
output=/your/output/directory

# Run featureCounts to quantify reads overlapping Batch1 peaks using expt2 BAM files
featureCounts -p -F SAF \
  -a ${saf_base}All_Samples.fwp.filter.non_overlapping.saf \
  --fracOverlap 0.2 \
  -o ${output}b1_on_b2_All_Samples_summit.counts \
  ${bam_base}Bulk_ATAC_21.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_22.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_23.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_24.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_25.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_26.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_27.noMT.filtered.deduped.bam \
  ${bam_base}Bulk_ATAC_28.noMT.filtered.deduped.bam