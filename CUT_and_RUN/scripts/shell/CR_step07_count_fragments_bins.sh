#!/bin/bash

###############################################################################
# Script: CR_step07_count_fragments_bins.sh
# Purpose: For each sample, tabulate the number of fragments whose midpoints fall
#          within specified genomic bins (e.g. 500 bp) for visualization/workflow.
#
# Usage: 
#   1. Edit projPath and sample settings as needed.
#   2. Run: bash CR_step07_count_fragments_bins.sh
#
###############################################################################

# Optional: load modules for HPC
# module purge
# module load samtools
# module load bedtools

projPath="/path/to/project"

start_sample=1
end_sample=8       # Adjust for your sample count
binLen=500         # Bin width in bp (e.g. 500bp bins)
samples=()

for i in $(seq $start_sample $end_sample); do
    samples+=("CR6_$i")
done

for sample in "${samples[@]}"; do
    input_file="$projPath/${sample}_bowtie2.fragments.bed"
    output_file="$projPath/${sample}_bowtie2.fragmentsCount.bin${binLen}.bed"
    
    if [[ ! -f "$input_file" ]]; then
        echo "Warning: Input file missing: $input_file. Skipping $sample."
        continue
    fi

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Counting fragments in $sample..."

    # Bin fragment midpoints and count occurrences per bin
    awk -v w="$binLen" '{print $1, int(((($2+$3)/2)/w))*w + w/2}' "$input_file" | \
        sort -k1,1V -k2,2n | \
        uniq -c | \
        awk -v OFS="\t" '{print $2, $3, $1}' | \
        sort -k1,1V -k2,2n > "$output_file"

    echo "$(date '+%Y-%m-%d %H:%M:%S'): Finished $sample. Bin counts in $output_file"
done

echo "Fragment bin counting complete for all samples."