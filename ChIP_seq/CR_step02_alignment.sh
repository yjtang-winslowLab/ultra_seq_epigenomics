#!/bin/bash

## Load required modules
module purge
module load bowtie2
module load samtools

# Define parameters
cores=24
ref="/your/directory/Reference_genomes/Bowtie2_mm10/mm10"
projPath="/your/directory/"
fqPath="${projPath}/keuh_trim"
alignPath="${projPath}/keuh_alignment"

# Define the list of samples to process
samples=(
    "SRR22028893"
    "SRR22028897"
    "SRR22028898"
    "SRR22028899"
    "SRR22028900"
)

# Create necessary output directories
mkdir -p "${alignPath}/bam"
mkdir -p "${alignPath}/logs"

# Alignment function
run_alignment() {
    local sample=$1

    echo "$(date): Starting alignment for ${sample}..."

    local input_fq="${fqPath}/${sample}.trim.fq.gz"
    local output_bam="${alignPath}/bam/${sample}.sorted.bam"
    local log_file="${alignPath}/logs/${sample}_bowtie2.log"

    # Check if input file exists
    if [[ ! -f "${input_fq}" ]]; then
        echo "Error: Input file ${input_fq} does not exist!" >&2
        return 1
    fi

    # Run Bowtie2 alignment, convert to BAM, and sort
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 \
        -I 10 -X 700 -p "${cores}" -x "${ref}" -U "${input_fq}" 2> "${log_file}" | \
        samtools view -@ "${cores}" -bS - | \
        samtools sort -@ "${cores}" -o "${output_bam}"

    # Check if alignment and sorting were successful
    if [[ $? -ne 0 ]]; then
        echo "Error: Alignment or sorting failed for sample ${sample}." >&2
        return 1
    fi

    # Index the sorted BAM file
    samtools index "${output_bam}"

    echo "$(date): Alignment, sorting, and indexing complete for ${sample}."
}

# Loop through samples and run alignment
for sample in "${samples[@]}"; do
    run_alignment "${sample}"
done

echo "$(date): All alignments completed successfully."