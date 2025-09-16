#!/bin/sh

## Replace MY_PI_SUNetID_or_Project_ID with the PI/project to be charged.
#SBATCH --account=mwinslow

## Set job time to 1 hour.
#SBATCH --time=1:00:00

## Send Email to me when done
#SBATCH --mail-user=yjtang@stanford.edu

## Set a name for the job, visible in `squeue`
#SBATCH --job-name=Peak_Annotation

## Quality of Service (QOS); think of it as sending your job into a special queue
#SBATCH --qos=normal

## One node.
#SBATCH --nodes=1

## One CPU/core per task
#SBATCH --cpus-per-task=4

## 128GB of RAM
#SBATCH --mem=128G

## Specify log file location to help with better logging of errors/outputs.

#SBATCH -o UA-%A-%a.out
#SBATCH -e UA-%A-%a.err

# Load required modules (adjust according to your HPC environment)
module purge
module load homer/4.11

# Annotating peaks with HOMER (annotatePeaks.pl)

# Define input directory (edit this path as needed)
input_base="path/to/your/input_directory/"

# Define output directory (edit this path as needed)
output_base="path/to/your/output_directory/"

# List of datasets to annotate (edit this array with your filenames, without extensions)
datasets=("b1_on_b2_peaks") 

# Loop through each dataset and run HOMER annotation
for dataset in "${datasets[@]}"
do
    echo "Running HOMER annotatePeaks on ${dataset}"

    # Define input and output filenames
    input_file="${input_base}${dataset}.txt"
    output_file="${output_base}${dataset}_annotated_peaks.txt"

    # Run HOMER annotatePeaks.pl
    annotatePeaks.pl "$input_file" mm10 > "$output_file"

    echo "Annotation complete for ${dataset}. Results saved to ${output_file}"
done