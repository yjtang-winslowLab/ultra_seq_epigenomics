#!/bin/sh

# Load necessary modules
module purge
module load samtools
module load bedtools

# Define paths
projPath="/your/directory/Keuh"
bamPath="${projPath}/keuh_alignment/bam"
bedPath="${projPath}/keuh_bed"
outputDir="${projPath}/keuh_bedgraph"
chromSize="/your/directory/Reference_genomes/mm10.chrom.sizes"
logPath="${projPath}/logs/bedgraph_logs"

# Create output and log directories
mkdir -p "$outputDir"
mkdir -p "$logPath"

# Effective genome size for mouse mm10
effectiveGenomeSize=2652783500

# Fragment length (average fragment size). Adjust if known from your data.
fragmentLength=100

# Define samples and labels
declare -A samples=(
    ["SRR22028893"]="igG"
    ["SRR22028897"]="H3K14ac_rep4"
    ["SRR22028898"]="H3K14ac_rep3"
    ["SRR22028899"]="H3K14ac_rep1"
    ["SRR22028900"]="H3K14ac_rep2"
)

# Loop through samples
for sample in "${!samples[@]}"; do
    sampleLabel="${samples[$sample]}"
    bamFile="${bamPath}/${sample}.sorted.bam"
    bedFile="${bedPath}/${sample}.bed"
    outputFile="${outputDir}/${sampleLabel}_bowtie2.fragments.norm.bedgraph"
    logFile="${logPath}/${sampleLabel}_bedgraph.log"

    {
        echo "Processing sample: ${sampleLabel} (${sample})"
        echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"

        # Check if BAM and BED files exist
        if [[ ! -f "$bamFile" ]]; then
            echo "Error: BAM file $bamFile not found!"
            continue
        fi
        if [[ ! -f "$bedFile" ]]; then
            echo "Error: BED file $bedFile not found!"
            continue
        fi

        # Count total mapped reads from BAM file
        totalReads=$(samtools view -c -F 260 "$bamFile")
        echo "Total mapped reads: $totalReads"

        # Check if totalReads is zero to avoid division by zero
        if [[ "$totalReads" -eq 0 ]]; then
            echo "Error: No mapped reads found in $bamFile!"
            continue
        fi

        # Calculate scaling factor
        scalingFactor=$(echo "scale=10; $effectiveGenomeSize / ($totalReads * $fragmentLength)" | bc)
        echo "Scaling factor: $scalingFactor"

        # Generate normalized bedgraph using bedtools genomecov
        echo "Generating normalized bedgraph..."
        bedtools genomecov -bg -scale "$scalingFactor" -i "$bedFile" -g "$chromSize" > "$outputFile"

        # Sort bedgraph file
        sort -k1,1 -k2,2n "$outputFile" -o "$outputFile"

        echo "Completed normalized bedgraph for ${sampleLabel}."
        echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
    } &> "$logFile"

done 

echo "$(date '+%Y-%m-%d %H:%M:%S'): All normalized bedgraph files generated."
