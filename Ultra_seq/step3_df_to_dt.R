################################################################################
# Ultra-Seq: Analysis Pipeline Step 3
# Convert master data frame to sample/gene/cell-number dictionary
#
# Author: Yuning J. Tang
# Affiliation: Dept of Genetics, Stanford University
#
# - Converts master data frame to a nested list (like a dictionary):
#   mouse/sample → gene/guide → cell numbers
################################################################################

library(tidyverse)

## --- FUNCTIONS --- ##

### Convert data frame to nested dictionary: sample → gene → cell numbers
df_to_cellnum_dict <- function(df) {
  split_by_sample <- split(df, df$sample_ID)
  # For each sample, split further by guide/gene
  dict <- lapply(
    split_by_sample, 
    function(sample_df) {
      split_by_gene <- split(sample_df, sample_df$guide)
      # For each gene, retain just the cell_num vector
      lapply(split_by_gene, function(gene_df) gene_df$cell_num)
    }
  )
  return(dict)
}

## --- RUN PIPELINE --- ##
cellnum_dict <- df_to_cellnum_dict(master_df_cel_final)

## --- SAVE OUTPUT --- ##
saveRDS(cellnum_dict, file = output_dictionary_rds)
message("Step 3 complete. Dictionary saved as: ", output_dictionary_rds)