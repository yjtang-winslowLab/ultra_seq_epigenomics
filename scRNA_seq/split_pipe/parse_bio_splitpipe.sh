# Purge modules to ensure a clean environment
module purge

# Activate Miniconda environment for Split-Pipe
source /your/directory/miniconda3/bin/activate

# (Optional: ensure conda is initialized in bash)
conda init bash

# Activate the specific conda environment containing Split-Pipe and dependencies
conda activate spipe

# Run Split-Pipe with merged FASTQ files and sample sheet
split-pipe \
    --mode all \
    --chemistry v2 \
    --genome_dir /labs/mwinslow/Reference_genomes/parse_mm10 \
    --fq1 /labs/mwinslow/Jackie/Parse_WT/merged_lib/merge_sub_lib1_R1.fastq.gz \
    --fq2 /labs/mwinslow/Jackie/Parse_WT/merged_lib/merge_sub_lib1_R2.fastq.gz \
    --output_dir /labs/mwinslow/Jackie/Parse_WT/output/2023_05_15_czi/sublib1 \
    --sample sgMeaf6_LA936 A1-A3 \
    --sample sgSafe14_MT2110 A4-A6 \
    --sample sgKmt2a_HC2555 A7-A9 \
    --sample sgPsip1_MT2140 A10-A12 \
    --sample sgZmynd8_LA498 B1-B3 \
    --sample sgPsip1_MT2131 B4-B6 \
    --sample sgKmt2a_HC2550 B7-B9 \
    --sample sgSafe14_MT2179 B10-B12 \
    --sample sgPsip1_MT2166 C1-C3 \
    --sample sgMeaf6_MT2176 C4-C6 \
    --sample sgSafe14_MT2144 C7-C9 \
    --sample sgZmynd8_MT1994 C10-C12 \
    --sample sgKmt2a_HC2551 D1-D3 \
    --sample sgMeaf6_MT2141 D4-D6 \
    --sample sgZmynd8_CM3581 D7-D9 \
    --sample sgSafe14_MT2179_2 D10-D12

# Deactivate the conda environment after processing
conda deactivate