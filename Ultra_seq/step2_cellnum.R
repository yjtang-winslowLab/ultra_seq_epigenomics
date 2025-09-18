################################################################################
# Ultra-Seq: Analysis Pipeline Step 2
# Generate cell number using spike-in
#
# Author: Yuning J. Tang
# Affiliation: Dept of Genetics, Stanford University
#
# 1. Calculates spike-in depth for each spike-in
# 2. Calculates cell number for each tumor based on spike-in depth
# 3. QC plots: spike-in depth, representation, number of tumors per sample
################################################################################

library(tidyverse)

spike_in_name <- c('spikein_1', 'spikein_2', 'spikein_3')

##### --- 1: Spike-in Depth Calculations --- #####
spike_depth <- master_df %>%
  filter(guide %in% spike_in_name) %>%
  group_by(sample_ID, guide) %>%
  summarise(spike_total = sum(bc_read_num), .groups="drop") %>%
  mutate(spike_read_per_cell = spike_total / 100000) # Note: In the sample prep step, we add 100,000 cells per spike-in cell line into each sample. 

spike_depth_mean <- spike_depth %>%
  group_by(sample_ID) %>%
  summarise(
    spike_depth_mean = mean(spike_read_per_cell),
    spike_depth_sd = sd(spike_read_per_cell),
    spike_depth_median = median(spike_read_per_cell),
    .groups="drop"
  )

spike_depth_mean_summary <- spike_depth %>%
  group_by(guide) %>%
  summarise(
    mean_depth = mean(spike_read_per_cell),
    sd_depth = sd(spike_read_per_cell),
    median_depth = median(spike_read_per_cell),
    .groups="drop"
  )

# --- QC PLOTS ---
spike_depth_plot <- ggplot(spike_depth, aes(x = reorder(factor(sample_ID), spike_read_per_cell), y = spike_read_per_cell, color = guide)) +
  geom_point() + ggtitle("Spike-in depth for each sample") +
  xlab("Sample ID") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(spike_depth_plot)

spike_depth_guide_plot <- ggplot(spike_depth, aes(x = guide, y = spike_read_per_cell, color = guide)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha = 0.8) +
  theme_classic() + ggtitle("Spike-in depth") +
  xlab("Spike-in vectors")
print(spike_depth_guide_plot)

spike_depth_comb_plot <- ggplot(spike_depth_mean, aes(x = 'spike_in depth', y = spike_depth_mean)) +
  geom_boxplot() +
  ggtitle("Spike-in depth (per sample mean)") +
  xlab("Spike-in vectors")
print(spike_depth_comb_plot)

##### --- 2: Calculate Cell Number Per Tumor --- #####
master_df_cel <- left_join(master_df, spike_depth_mean, by="sample_ID") %>%
  mutate(cell_num = bc_read_num / spike_depth_mean)

##### --- 3: Spike-in Representation --- #####
spike_mean <- master_df %>%
  filter(guide %in% spike_in_name) %>%
  group_by(sample_ID) %>%
  summarise(spike_total = sum(bc_read_num), .groups="drop") %>%
  mutate(spike_per_cell = spike_total / 300000)

spike_representation <- master_df %>%
  group_by(sample_ID) %>%
  summarise(total = sum(bc_read_num), .groups="drop") %>%
  left_join(spike_mean, by = "sample_ID") %>%
  mutate(spike_pct = spike_total / total)

spike_representation_plot <- ggplot(spike_representation, aes(x = reorder(factor(sample_ID), spike_pct), y = spike_pct, fill = sample_ID)) +
  geom_bar(stat="identity", width=0.5, alpha = 0.9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(legend.position = "none")
print(spike_representation_plot)

##### --- 4: QC: Tumor Number per Sample --- #####
tn <- master_df_cel %>%
  filter(cell_num >= 300, !guide %in% spike_in_name) %>%
  group_by(sample_ID) %>%
  summarise(n = n(), .groups="drop") %>%
  left_join(expt_info, by = "sample_ID")

tumor_count_plot <- ggplot(tn, aes(x = reorder(sample_ID, -n), y = n)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab('Total tumor number') +
  xlab('Sample ID')
print(tumor_count_plot)

##### --- 5: Remove Spike-ins & Output --- #####
# Remove spike-ins and specific guides from final dataframe
rm_guides <- c(spike_in_name, 'Cdkn1a_2', 'Setd2_2', 'Arid1a_2')
master_df_cel_final <- master_df_cel %>%
  filter(!guide %in% rm_guides)

message("Step 2 complete. Final cell-number dataframe written to: ", output_cel_final)
