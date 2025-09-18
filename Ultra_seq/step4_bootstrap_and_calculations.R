################################################################################
# Ultra-Seq: Analysis Pipeline Step 4
# Bootstrap and Calculate Tumor Metrics
#
# Author: Yuning J. Tang
# Affiliation: Dept of Genetics, Stanford University
#
# Functions:
#  - bootstrap_mouse
#  - bootstrap_tumors
#  - data_Prep (combine tumors per guide)
#  - ln_mean, calc_LNMean
#  - calc_percentile
#  - mouse_group (split KT/KTC for normalization)
#  - calc_TumorNum, calc_TumorBurden
################################################################################

# Dependencies: assumes tidyverse loaded
# User must supply:
#   guide_list: vector of guide names/genes
#   inert_id: vector of inert/control guide names
#   kt_ID: vector of KT mouse IDs
#   sample_num: number of mice/samples
#   l: number of guides in guide_list

#########################################################################
# -- 1. Bootstrap mouse selection (sample with replacement, enforce min KT)
#########################################################################
bootstrap_mouse <- function(simp_dict, kt_ID, required_num_KT) {
  repeat {
    bs_mouse_dict <- sample(simp_dict, length(simp_dict), replace = TRUE)
    kt_check <- sum(names(bs_mouse_dict) %in% kt_ID)
    if (kt_check >= required_num_KT) break
  }
  return(bs_mouse_dict)
}

#########################################################################
# -- 2. Bootstrap tumors within each mouse/sample
# Returns list: mouse > guide > cell number
#########################################################################
bootstrap_tumors <- function(bs_mouse_dict, guide_list, min_tumor_size) {
  sample_num <- length(bs_mouse_dict)
  boot_data <- vector("list", sample_num)
  names(boot_data) <- names(bs_mouse_dict)
  for (m in seq_len(sample_num)) {
    sample_id <- names(bs_mouse_dict)[m]
    boot_data[[m]] <- vector("list", length(guide_list))
    names(boot_data[[m]]) <- guide_list
    for (g in guide_list) {
      guide_data <- bs_mouse_dict[[m]][[g]]
      if (is.null(guide_data) || length(guide_data) == 0) {
        boot_data[[m]][[g]] <- NA
      } else {
        booted <- sample(guide_data, length(guide_data), replace = TRUE)
        # Keep only tumors above size threshold
        boot_data[[m]][[g]] <- booted[booted >= min_tumor_size]
      }
    }
  }
  return(boot_data)
}

#########################################################################
# -- 3. Combine all tumors from mice for each guide
# Returns list: guide > cell number vector (all samples combined)
#########################################################################
data_Prep <- function(boot_tumors, guide_list) {
  guide_data <- setNames(vector("list", length(guide_list)), guide_list)
  for (m in seq_along(boot_tumors)) {
    for (g in guide_list) {
      data <- boot_tumors[[m]][[g]]
      if (!is.null(data) && length(data) > 0) {
        guide_data[[g]] <- c(guide_data[[g]], data)
      }
    }
  }
  return(guide_data)
}

#########################################################################
# -- 4. Metrics Functions -------------------------------------------------
#########################################################################
# Log-normal mean
ln_mean <- function(cell_num) {
  ln_cm <- log(cell_num)
  x <- mean(ln_cm)
  vx <- var(ln_cm)
  exp(x + 0.5 * vx)
}

# Log-normal mean for all guides, normalized to inert
# Returns data.frame: guide, lnMean, Inert_lnMean, relative_ln_mean
calc_LNMean <- function(prepared_data, guide_list, inert_id) {
  res <- data.frame(guide = guide_list, lnMean = NA, Inert_lnMean = NA, relative_ln_mean = NA)
  for (g in seq_along(guide_list)) {
    cells <- na.omit(prepared_data[[g]])
    res$lnMean[g] <- if (length(cells) > 1) ln_mean(cells) else NA
  }
  # Inert calculation
  inert_means <- res$lnMean[res$guide %in% inert_id]
  inert_median <- median(inert_means, na.rm = TRUE)
  res$Inert_lnMean <- inert_median
  res$relative_ln_mean <- res$lnMean / inert_median
  return(res)
}

# Percentile sizes for each guide (normalized to inert median for each percentile)
# pct = vector of percentiles (e.g. c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95))
calc_percentile <- function(prepared_data, guide_list, inert_id, pct = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)) {
  res <- data.frame(guide = guide_list)
  for (g in seq_along(guide_list)) {
    cells <- na.omit(prepared_data[[g]])
    res[g, 2:(length(pct)+1)] <- if (length(cells) > 1) quantile(cells, pct) else rep(NA, length(pct))
  }
  colnames(res)[2:(length(pct)+1)] <- paste0("pct_size_", pct)
  # Normalize by inert median for each percentile
  inert_rows <- res$guide %in% inert_id
  inert_medians <- apply(res[inert_rows,-1, drop=FALSE], 2, median, na.rm = TRUE)
  norm_data <- sweep(res[,-1, drop=FALSE], 2, inert_medians, "/")
  res[,-1] <- norm_data
  return(res)
}

#########################################################################
# -- 5. Split out KT vs. KTC mouse/sample groups
#########################################################################
# Returns list: ktc_tumors and kt_tumors (each guide > cell numbers)
mouse_group <- function(boot_tumors, guide_list, kt_ID) {
  kt_list <- list()
  ktc_list <- list()
  for (m in seq_along(boot_tumors)) {
    sample_id <- names(boot_tumors)[m]
    if (sample_id %in% kt_ID) {
      kt_list[[length(kt_list)+1]] <- boot_tumors[[m]]
    } else {
      ktc_list[[length(ktc_list)+1]] <- boot_tumors[[m]]
    }
  }
  # Combine across mice:
  out <- list(
    ktc_tumors = data_Prep(ktc_list, guide_list),
    kt_tumors = data_Prep(kt_list, guide_list)
  )
  return(out)
}

#########################################################################
# -- 6. Tumor Number: calculate for each guide, normalized
#########################################################################
calc_TumorNum <- function(grouped_list, guide_list, inert_id) {
  ktc <- grouped_list[["ktc_tumors"]]
  kt <- grouped_list[["kt_tumors"]]
  
  df <- data.frame(guide = guide_list,
                   ktc_number_raw = sapply(ktc, function(x) sum(!is.na(x))),
                   kt_number_raw = sapply(kt, function(x) sum(!is.na(x))),
                   tumor_number_norm = NA,
                   relative_tumor_num = NA)
  df$tumor_number_norm <- df$ktc_number_raw / df$kt_number_raw
  inert_norms <- df$tumor_number_norm[df$guide %in% inert_id]
  inert_median <- median(inert_norms, na.rm = TRUE)
  df$relative_tumor_num <- df$tumor_number_norm / inert_median
  return(df)
}

#########################################################################
# -- 7. Tumor Burden: calculate for each guide, normalized
#########################################################################
calc_TumorBurden <- function(grouped_list, guide_list, inert_id) {
  ktc <- grouped_list[["ktc_tumors"]]
  kt <- grouped_list[["kt_tumors"]]
  
  df <- data.frame(guide = guide_list,
                   ktc_burden_raw = sapply(ktc, function(x) sum(x, na.rm = TRUE)),
                   kt_burden_raw = sapply(kt, function(x) sum(x, na.rm = TRUE)),
                   tumor_burden_norm = NA,
                   relative_tumor_burden = NA)
  df$tumor_burden_norm <- df$ktc_burden_raw / df$kt_burden_raw
  inert_norms <- df$tumor_burden_norm[df$guide %in% inert_id]
  inert_median <- median(inert_norms, na.rm = TRUE)
  df$relative_tumor_burden <- df$tumor_burden_norm / inert_median
  return(df)
}

################################################################################
# End of script. All functions above are ready for modular use and sharing.
# Typical workflow:
# bs_mouse <- bootstrap_mouse(simp_dict, kt_ID)
# boot_tumors <- bootstrap_tumors(bs_mouse, guide_list, min_tumor_size)
# guide_data <- data_Prep(boot_tumors, guide_list)
# ln_means <- calc_LNMean(guide_data, guide_list, inert_id)
# pct_sizes <- calc_percentile(guide_data, guide_list, inert_id)
# mouse_split <- mouse_group(boot_tumors, guide_list, kt_ID)
# tumor_nums <- calc_TumorNum(mouse_split, guide_list, inert_id)
# tumor_burdens <- calc_TumorBurden(mouse_split, guide_list, inert_id)
################################################################################