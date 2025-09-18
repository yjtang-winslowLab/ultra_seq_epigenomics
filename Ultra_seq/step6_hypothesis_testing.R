################################################################################
# Ultra-Seq: Step 6 - Compute 95% confidence intervals, p-values, FDR for metrics,
#                   gene-level scores, and combine p-values by Stouffer's method.
#
# Author: [Your Name]
# Date: [YYYY-MM-DD]
################################################################################

library(tidyverse)

# ===== CONFIGURATION =====
# Update these to match your environment
output_folder <- "output/tables/point_est/KTC/"
kt_ID <- c('S9','S15','S30')

metric_names <- c("log_normal_mean",
                  "tumor_number",
                  "tumor_burden",
                  "percentile_size_0.5", "percentile_size_0.6",
                  "percentile_size_0.7", "percentile_size_0.8",
                  "percentile_size_0.9", "percentile_size_0.95")

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

dict_for_stats <- metrics_bootstrap
guide_list <- names(dict_for_stats)
l <- length(guide_list)

# ===== FUNCTIONS =====

fill_data <- function(stat_dict, metric){
  data.frame(
    guide    = guide_list,
    point_est= sapply(guide_list, function(g) mean(as.numeric(stat_dict[[g]][[metric]]), na.rm=TRUE))
  )
}

ci_95 <- function(stat_dict, metric) {
  ci <- c(0.025, 0.975)
  data.frame(
    guide    = guide_list,
    lower_ci = sapply(guide_list, function(g) quantile(na.omit(as.numeric(stat_dict[[g]][[metric]])), ci, na.rm=TRUE)[1]),
    upper_ci = sapply(guide_list, function(g) quantile(na.omit(as.numeric(stat_dict[[g]][[metric]])), ci, na.rm=TRUE)[2])
  )
}

sig_test <- function(stat_dict, metric){
  mu0 <- 1
  data.frame(
    guide    = guide_list,
    p_onetail= sapply(guide_list, function(g) {
      nums <- as.numeric(stat_dict[[g]][[metric]])
      min(sum(nums <= mu0), sum(nums > mu0)) / length(nums)
    }),
    p_twotail= sapply(guide_list, function(g) {
      nums <- as.numeric(stat_dict[[g]][[metric]])
      min_alpha <- min(sum(nums <= mu0)/length(nums), sum(nums > mu0)/length(nums))
      min_alpha * 2
    })
  ) %>% mutate(fdr = p.adjust(p_twotail, method = 'fdr'))
}

assign_gene <- function(df){
  gene_symbol <- sapply(df$guide, function(id) {
    if (substr(id, 1,2) == "sg") id else str_to_title(str_extract(id, "^.+?(?=_)"))
  })
  add_column(df, gene_name = gene_symbol, .after = 2)
}

#################################################################
##### weighted z method for gene level metrics ###################
# ref: Chen, PMID: 21401770; Zaykin, PMID: 21605215; Whitlock, PMID: 16135132
# Steps:
# 1. convert p value to z score
# 2. assign weights to each z score. 
#    # Weights in Stouffer's method are assigned as the number of tumors observed in the KTC genotype > 300 cells 
#    # (each sgRNA's weight is its tumor count among KTC mice)
stouffer.test <- function(p_onetail, w) { 
  p_onetail[p_onetail == 1] <- 0.99999
  z_score <- -qnorm(p_onetail)
  z_weighted <- sum(w*z_score)/sqrt(sum(w^2))
  stouffer_alpha <- 1 - pnorm(z_weighted)
  stouffer_alpha_min <- min(stouffer_alpha, 1 - stouffer_alpha)
  stouffer_p <- stouffer_alpha_min * 2
  return(stouffer_p)
}
#################################################################

gene_score <- function(df, w) {
  df$gene_effect <- df$point_est * w$pct_weight
  df %>% group_by(gene_name) %>% summarise(point_est = sum(gene_effect), .groups = "drop")
}

gene_stat <- function(df, w){
  df$z_weight <- w$tumor_number
  stouffer_p <- df %>%
    group_by(gene_name) %>%
    summarise(p = stouffer.test(p_onetail, z_weight), .groups = "drop")
  stouffer_p$fdr <- p.adjust(stouffer_p$p, method = "fdr", nrow(stouffer_p))
  outs <- gene_score(df, w)
  left_join(outs, stouffer_p, by = "gene_name")
}

# ===== CALCULATIONS =====

# 1. Metrics per guide
metrics_raw <- lapply(metric_names, function(m) fill_data(dict_for_stats, m))
names(metrics_raw) <- metric_names

ci_list <- lapply(metric_names, function(m) ci_95(dict_for_stats, m))
names(ci_list) <- metric_names

pval_list <- lapply(metric_names, function(m) sig_test(dict_for_stats, m))
names(pval_list) <- metric_names

metrics_df_list <- list()
for (m in metric_names) {
  metrics_df_list[[m]] <- assign_gene(
    bind_cols(metrics_raw[[m]], 
              ci_list[[m]][, c("lower_ci", "upper_ci")], 
              pval_list[[m]][, c("p_onetail", "p_twotail", "fdr")])
  )
}

# 2. Weights (include all guides)
# Get complete list of all guides from the dataset:
all_guides <- unique(master_df_cel_final$guide)

weights_raw <- master_df_cel_final %>%
  select(guide, cell_num, Geotype) %>%
  filter(cell_num >= 300 & Geotype == "KTC") %>%
  group_by(guide) %>%
  summarise(tumor_number = n(), .groups = "drop")

# Add zero entries for guides with no qualifying tumors, so ALL guides are present:
weights_complete <- tibble(guide = all_guides) %>%
  left_join(weights_raw, by = "guide") %>%
  mutate(tumor_number = replace_na(tumor_number, 0))

# Compute pct_weight for each gene
weights_complete <- assign_gene(weights_complete)
gene_tumor_count <- weights_complete %>%
  group_by(gene_name) %>%
  summarise(total = sum(tumor_number), .groups = "drop")
weights_complete <- left_join(weights_complete, gene_tumor_count, by = "gene_name")
weights_complete$pct_weight <- ifelse(weights_complete$total > 0, weights_complete$tumor_number / weights_complete$total, 0)

# 3. Gene-level summary using the complete weights
gene_score_ln <- gene_stat(metrics_df_list[["log_normal_mean"]], weights_complete)
gene_score_tn <- gene_stat(metrics_df_list[["tumor_number"]], weights_complete)
gene_score_tb <- gene_stat(metrics_df_list[["tumor_burden"]], weights_complete)
gene_score_0.9 <- gene_stat(metrics_df_list[["percentile_size_0.9"]],weights_complete)

# 4. Export tables -- gene-level
write.csv(gene_score_ln,  file.path(output_folder, "gene_score_ln_KTC.csv"), row.names=FALSE)
write.csv(gene_score_tn,  file.path(output_folder, "gene_score_tn_KTC.csv"), row.names=FALSE)
write.csv(gene_score_tb,  file.path(output_folder, "gene_score_tb_KTC.csv"), row.names=FALSE)
write.csv(gene_score_0.9, file.path(output_folder, "gene_score_0.9_KTC.csv"), row.names=FALSE)

# 5. Export sgRNA-level tables -- all metrics
for (m in names(metrics_df_list)) {
  write.csv(metrics_df_list[[m]], file.path(output_folder, paste0(m, "_KTC.csv")), row.names=FALSE)
}

# 6. All percentiles combined (long format)
pct_metrics_names <- metric_names[grepl("^percentile_size_", metric_names)]
pct_final_df <- bind_rows(lapply(pct_metrics_names, function(m) metrics_df_list[[m]]))
write.csv(pct_final_df, file.path(output_folder, "pct_final_df_KTC.csv"), row.names=FALSE)

message("Step 6 complete. Results saved to ", output_folder)