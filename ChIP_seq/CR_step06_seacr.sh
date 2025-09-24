#!/bin/sh

# Load necessary modules
module purge
module load SEACR/1.3
module load bedtools
module load R

# Define project path for input files
projPath="/your/directory/"
seacr="/your/apps/software/SEACR/1.3/SEACR_1.3.sh"

# Define samples
files=(
    "H3K14ac_rep1"
    "H3K14ac_rep2"
    "H3K14ac_rep3"
    "H3K14ac_rep4"
)

# Define IgG control sample
controlSample="igG"

# Create necessary directories for SEACR output and logs
outputDir="${projPath}/keuh_SEACR"
logPath="${projPath}/logs/seacr_logs"

mkdir -p "${outputDir}"
mkdir -p "${logPath}"

# Check if SEACR script exists and is executable
if [[ ! -x "${seacr}" ]]; then
    echo "Error: SEACR script ${seacr} not found or not executable!" >&2
    exit 1
fi

# Loop through each target sample for SEACR peak calling
for file in "${files[@]}"; do

    inputBed="${projPath}/keuh_bedgraph/${file}_bowtie2.fragments.norm.bedgraph"
    controlBed="${projPath}/keuh_bedgraph/${controlSample}_bowtie2.fragments.norm.bedgraph"
    logFile="${logPath}/${file}_seacr.log"

    {
        echo "---------------------------------------------"
        echo "$(date): Starting SEACR peak calling for sample ${file}"

        # Check if input BED files exist
        if [[ ! -f "${inputBed}" ]]; then
            echo "Error: Input BED file ${inputBed} not found!" >&2
            continue
        fi

        if [[ ! -f "${controlBed}" ]]; then
            echo "Error: Control BED file ${controlBed} not found!" >&2
            continue
        fi

        # Run SEACR for histone peaks against IgG control
        echo "$(date): Running SEACR..."
        bash "${seacr}" "${inputBed}" "${controlBed}" norm stringent "${outputDir}/${file}_seacr_control.peaks"

        if [[ $? -ne 0 ]]; then
            echo "Error: SEACR failed for sample ${file}!" >&2
            continue
        fi

        echo "$(date): SEACR peak calling completed successfully for sample ${file}"

    } &> "${logFile}"

done

echo "$(date): SEACR peak calling completed for all samples."
