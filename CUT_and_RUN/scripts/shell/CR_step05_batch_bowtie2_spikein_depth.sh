#!/bin/bash

###############################################################################
# Script: CR_step05_batch_bowtie2_spikein_depth.sh
# Purpose: Align trimmed paired-end reads to E. coli spike-in reference using Bowtie2,
#          compute sequencing depth (mapped reads), and save results per sample.
#
# Usage:
#   1. Edit ref, projPath, fqPath, core count, sample range below.
#   2. Run with: bash CR_step05_batch_bowtie2_spikein_depth.sh
###############################################################################

# For cluster environments: load Bowtie2 & Samtools modules
# module purge
# module load bowtie2/2.2.0
# module load samtools

set -e  # Stop on first error

cores=24
ref="/path/to/Reference_genomes/ecoli/ecoli"          # Bowtie2 index basename for E. coli spike-in
projPath="/path/to/project"                           # Project root/output base
fqPath="/path/to/trimmed_reads"                       # Directory with trimmed reads

# Sample naming and count
start_sample=1
end_sample=8     # Change if you have more/fewer samples

samples=()
for i in $(seq "$start_sample" "$end_sample"); do
    samples+=("CR_$i")
done

# Output locations
samDir="${projPath}/alignment/ecoli_sam"
summaryDir="${samDir}/bowtie2_Spikein_summary"
mkdir -p "$samDir" "$summaryDir"

# Loop for Bowtie2 mapping and depth stats
for sample in "${samples[@]}"; do
    echo "$(date): Aligning $sample to E. coli (spike-in)..."

    fq1="${fqPath}/merged_${sample}_1.paired.trim.fq.gz"
    fq2="${fqPath}/merged_${sample}_2.paired.trim.fq.gz"
    samOut="$samDir/${sample}_spikeIn.sam"
    summaryOut="$summaryDir/${sample}_bowtie2_spikeIn.txt"

    # Bowtie2 alignment
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p "$cores" \
        -x "$ref" -1 "$fq1" -2 "$fq2" -S "$samOut" \
        &> "$summaryOut"

    # Check if aligner ran successfully
    if [[ $? -ne 0 ]]; then
        echo "Error occurred during Bowtie2 alignment for sample $sample"
        continue
    fi

    # Confirm SAM file existence
    if [[ ! -f "$samOut" ]]; then
        echo "Error: SAM file for $sample not found!" >&2
        continue
    fi

    # Calculate mapped read count (paired-end; samtools counts alignments, so divide by 2)
    echo "Calculating sequencing depth for $sample..."
    seqDepthDouble=$(samtools view -c -F 4 "$samOut")   # Only mapped reads

    seqDepth=$((seqDepthDouble / 2))
    echo "Total mapped reads for $sample: $seqDepthDouble"
    echo "Effective sequencing depth (pairs): $seqDepth"

    # Save depth value to file
    echo "$seqDepth" > "$summaryDir/${sample}_bowtie2_spikeIn.seqDepth"

    echo "Bowtie2 alignment and depth calculation for $sample complete."
done

echo "All alignments and depth calculations complete for samples $start_sample through $end_sample."