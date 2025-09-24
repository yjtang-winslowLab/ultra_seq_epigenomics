#!/bin/bash

module purge
module load samtools
module load bedtools

projPath="/your/directory/Keuh"
bamPath="${projPath}/keuh_alignment/bam"
bedPath="${projPath}/keuh_bed"
logPath="${projPath}/logs/bam_bed_logs"

mkdir -p "$bedPath"
mkdir -p "$logPath"

cores=8

samples=(
    "SRR22028893"
    "SRR22028897"
    "SRR22028898"
    "SRR22028899"
    "SRR22028900"
)

for sample_name in "${samples[@]}"; do
    sorted_bam="${bamPath}/${sample_name}.sorted.bam"
    bed_file="${bedPath}/${sample_name}.bed"
    logFile="${logPath}/${sample_name}_bed.log"

    {
        echo "---------------------------------------------"
        echo "$(date): Starting processing for sample ${sample_name}"

        if [[ ! -f "$sorted_bam" ]]; then
            echo "Error: Sorted BAM file $sorted_bam not found!" >&2
            continue
        fi

        # Optional: Index BAM if needed downstream
        if [[ ! -f "${sorted_bam}.bai" ]]; then
            echo "$(date): Indexing BAM file ${sorted_bam}..."
            samtools index -@ "$cores" "$sorted_bam"
            if [[ $? -ne 0 ]]; then
                echo "Error: Indexing failed for ${sorted_bam}!" >&2
                continue
            fi
        else
            echo "$(date): BAM index already exists for ${sorted_bam}."
        fi

        # Convert sorted BAM to BED format (single-end)
        echo "$(date): Converting BAM to BED for ${sample_name}..."
        bedtools bamtobed -i "$sorted_bam" | sort -k1,1 -k2,2n > "$bed_file"

        if [[ $? -ne 0 ]]; then
            echo "Error: BED conversion failed for ${sample_name}!" >&2
            continue
        fi

        echo "$(date): Successfully created BED file for ${sample_name}."

    } &> "$logFile"
done

echo "$(date): All samples processed successfully."