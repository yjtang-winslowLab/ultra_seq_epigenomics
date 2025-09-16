#!/bin/sh

# Load required modules (adjust according to your HPC environment)
module purge
module load homer/4.11

## Finding Enriched Motifs in Genomic Regions using HOMER (findMotifsGenome.pl)

## Define input directory (edit this path as needed)
input_base="path/to/your/input_directory"

## Define output directory (edit this path as needed)
output_base="path/to/your/ouput_directory"

## List of datasets (edit this array with your filenames, without extensions)
datasets=(kat7_bed_dn kat7_bed_up)

## Path to HOMER genome data (edit this path as needed)
genome_path="path/to/your/genome_directory"

## Run HOMER findMotifsGenome.pl for each dataset
for dataset in "${datasets[@]}"
do
    echo "Running HOMER findMotifsGenome on ${dataset}"

    # Create output directory for each dataset
    mkdir -p "${output_base}${dataset}"

    # Run HOMER motif enrichment analysis
    findMotifsGenome.pl "${input_base}${dataset}.txt" "${genome_path}" "${output_base}${dataset}" -size 200

    echo "Motif analysis complete for ${dataset}. Results saved to ${output_base}${dataset}"
done