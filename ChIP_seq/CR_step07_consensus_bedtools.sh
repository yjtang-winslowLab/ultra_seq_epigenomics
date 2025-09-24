#!/bin/sh

## Replace MY_PI_SUNetID_or_Project_ID with the PI/project to be charged.
#SBATCH --account=mwinslow

## Set job time to 1 hour.
#SBATCH --time=0:30:00

## Send Email to me when done
#SBATCH --mail-user=yjtang@stanford.edu

## Set a name for the job, visible in `squeue`
#SBATCH --job-name=consensus

## Quality of Service (QOS); think of it as sending your job into a special queue
#SBATCH --qos=normal

## One node.
#SBATCH --nodes=1

## One CPU/core per task
#SBATCH --cpus-per-task=4

## RAM
#SBATCH --mem=16G

## Specify log file location to help with better logging of errors/outputs.

#SBATCH -o UA-%A-%a.out
#SBATCH -e UA-%A-%a.err

# Load the bedtools module
module purge
module load bedtools

# Define project path
projPath="/labs/mwinslow/Jackie/CUT_RUN/Keuh/"
outputPath="$projPath/keuh_consensus"

# Create necessary output directories if they do not exist
mkdir -p "$outputPath"

# Define the input BED files for H3K14ac
inputFiles=(
    "$projPath/keuh_SEACR/H3K14ac_rep1_seacr_control.peaks.stringent.bed"
    "$projPath/keuh_SEACR/H3K14ac_rep2_seacr_control.peaks.stringent.bed"
    "$projPath/keuh_SEACR/H3K14ac_rep3_seacr_control.peaks.stringent.bed"
    "$projPath/keuh_SEACR/H3K14ac_rep4_seacr_control.peaks.stringent.bed"
)

# Perform the intersection for H3K14ac using bedtools 
bedtools intersect -a "${inputFiles[0]}" -b "${inputFiles[1]}" "${inputFiles[2]}" "${inputFiles[3]}" \
> "$outputPath/consensus_keuh_H3K14ac.bed"

echo "Consensus regions for H3K14ac written to $outputPath/consensus_keuh_H3K14ac.bed"
