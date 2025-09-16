#!/bin/bash

# Load required modules (adjust according to your HPC environment)
module purge
module load bedops/2.4.41

# Define base directory (edit this path as needed)
base="path/to/your/output_directory"

## Format the merged BED files into SAF format for downstream analysis
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="all_samples_summits_"++nr;  print peakid,$1,$2,$3,"."}' ${base}All_Samples.fwp.filter.non_overlapping.bed > ${base}All_Samples.fwp.filter.non_overlapping.saf

## Format the merged BED files with peak IDs for downstream analysis
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="all_samples_summits_"++nr;  print $1,$2,$3,peakid,"0","."}' ${base}All_Samples.fwp.filter.non_overlapping.bed > ${base}All_Samples.fwp.filter.non_overlapping.peakIds.bed