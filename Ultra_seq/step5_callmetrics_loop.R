################################################################################
# Ultra-Seq: Step 5 - Serial Bootstrap Metric Calculation - Single Loop
#
# Author: Yuning J. Tang
# Affiliation: Dept of Genetics, Stanford University
#
# This script:
#   - Runs two-step bootstrap sampling
#   - Computes metrics: log-normal mean, percentiles, tumor number, tumor burden
#   - Processes all guides over many bootstrap iterations
#   - Designed with generic KT mouse IDs for easy adaptation
################################################################################

library(tidyverse)

### ------------ CONFIGURATION ----------------

num_bootstrap <- 10
min_size <- 300
percentile <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

# ---- Specify KT genotype sample IDs ----
kt_ID <- c('S9','S15','S30')
required_num_KT <- length(kt_ID)

# -------- LOAD DATA AND PREPARE VARIABLES --------------

guide_list <- sort(unique(master_df_cel_final$guide))
l <- length(guide_list)
inert_id <- guide_list[grep('sg', guide_list)]
spike_ins <- guide_list[grep('spike', guide_list)]

metric_names <- c('log_normal_mean', 'tumor_number', 'tumor_burden',
                  paste0('percentile_size_', percentile))

# Pre-allocate output structure
metrics_bootstrap <- setNames(
  lapply(guide_list, function(g) {
    setNames(lapply(metric_names, function(m) numeric(num_bootstrap)), metric_names)
  }),
  guide_list
)

### -------------- MAIN SERIAL BOOTSTRAP LOOP ---------------

for (itr in seq_len(num_bootstrap)) {
  if (itr %% 10 == 0) {
    message("Bootstrap iteration: ", itr, "/", num_bootstrap)
  }
  
  boot_m <- bootstrap_mouse(cellnum_dict, kt_ID, required_num_KT)
  boot_t <- bootstrap_tumors(boot_m, guide_list, min_size)
  prepped_data <- data_Prep(boot_t, guide_list)
  groups <- mouse_group(boot_t, guide_list, kt_ID)
  # Compute metrics for KTC mice (change `$ktc_tumors` to `$kt_tumors` if desired)
  log_normal_mean <- calc_LNMean(groups$ktc_tumors, guide_list, inert_id)
  percentile_size <- calc_percentile(groups$ktc_tumors, guide_list, inert_id, percentile)
  tumor_number <- calc_TumorNum(groups, guide_list, inert_id)
  tumor_burden <- calc_TumorBurden(groups, guide_list, inert_id)
  
  for (g in seq_len(l)) {
    guide <- guide_list[g]
    metrics_bootstrap[[guide]][['log_normal_mean']][itr] <- log_normal_mean$relative_ln_mean[g]
    metrics_bootstrap[[guide]][['tumor_number']][itr] <- tumor_number$relative_tumor_num[g]
    metrics_bootstrap[[guide]][['tumor_burden']][itr] <- tumor_burden$relative_tumor_burden[g]
    for (p in seq_along(percentile)) {
      met_name <- paste0('percentile_size_', percentile[p])
      metrics_bootstrap[[guide]][[met_name]][itr] <- percentile_size[g, 1 + p]
    }
  }
}

### -------------- SAVE OUTPUT -----------------
saveRDS(metrics_bootstrap, "metrics_bootstrap_step5_serial.rds")
message("Step 5 complete. Bootstrapped metrics saved to: metrics_bootstrap_step5_serial.rds")