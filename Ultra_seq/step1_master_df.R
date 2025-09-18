################################################################################
# Ultra-Seq: Massively Parallel Quantification of Growth Phenotype by Sequencing
# Step 1: Create Master Data Frame
#
# Author: Yuning J. Tang
# Affiliation: Dept of Genetics, Stanford University
#
# This script reads processed data from multiple samples (.csv/.csv.gz),
# combines them, adds experimental info, and writes out a master data frame.
################################################################################

library(tidyverse)

# ---- Configuration (Update as needed!) ----
data_dir <- "/your/data/directory"   # Relative to your script location
expt_info_file <- "expt_info.csv"
output_file <- "output/master_df_step1.csv"

# ---- Read all sample CSV and CSV.GZ files in data_dir ----
scan_dir <- file.path(data_dir, "pilot_data")
sample_files <- list.files(scan_dir, pattern = "\\.csv$|\\.csv\\.gz$", full.names = TRUE)
if (length(sample_files) == 0) stop(paste("No CSV or CSV.GZ files found in", scan_dir))

# Read all files into list of data frames
data_list <- lapply(sample_files, read.csv)

# Extract sample IDs from filenames (before first underscore or dot)
sample_ids <- basename(sample_files) %>%
  str_replace("_.*", "") %>%
  str_replace("\\.csv(\\.gz)?$", "")

names(data_list) <- sample_ids

master_df <- bind_rows(data_list, .id = "sample_ID")

# ---- Experimental info merge ----
expt_info <- read.csv(file.path(data_dir, expt_info_file))
if (!"sample_ID" %in% colnames(expt_info)) {
  stop("Experimental info file must have a 'sample_ID' column.")
}
master_df <- left_join(master_df, expt_info, by = "sample_ID")

# ---- Write output ----
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_csv(master_df, output_file)

message("Step 1 complete. Master data frame written to: ", output_file)