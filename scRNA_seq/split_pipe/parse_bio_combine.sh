# Purge modules to ensure a clean environment
module purge

# Activate your Miniconda installation
source /path/to/miniconda3/bin/activate

# (Optional: run once per session)
conda init bash

# Activate the Conda environment containing Split-Pipe
conda activate spipe

# Define the base output directory for sublibraries
dir=/path/to/sublibrary/outputs/

# Combine outputs from multiple sublibraries into a single directory
split-pipe \
    --mode comb \
    --sublibraries ${dir}sublib1 ${dir}sublib2 ${dir}sublib4 ${dir}sublib5 ${dir}sublib6 ${dir}sublib7 ${dir}sublib8 \
    --output_dir /path/to/combined_output/no_lib3

# Deactivate the conda environment when finished
conda deactivate