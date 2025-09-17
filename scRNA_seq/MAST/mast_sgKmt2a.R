#!/usr/bin/env Rscript

# ============================================================
# MAST Single-cell Differential Expression Analysis (Parallel, Improved)
# Author: Yuning J. Tang
# Affiliation: Winslow Lab, Dept of Genetics, Stanford University
# ============================================================

# --- Load required libraries ---
message("Loading libraries...")

cran_pkgs <- c("tidyverse", "data.table")
bioc_pkgs <- c("SingleCellExperiment", "MAST", "BiocParallel")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

invisible(lapply(c(cran_pkgs, bioc_pkgs), function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# --- Set random seed for reproducibility ---
set.seed(12345)

# --- Set directories ---
data_dir <- "your/directory/data/"
output_dir <- "your/directory/output/"
if (!dir.exists(data_dir)) stop("Data directory does not exist: ", data_dir)
if (!dir.exists(output_dir)) dir.create(output_dir)
message("Using data dir: ", data_dir)
message("Using output dir: ", output_dir)

# --- Read metadata & counts ---
message("Reading cell metadata ...")
meta <- read.csv(file.path(data_dir, "kmt2a_filtered_metadata.csv"), stringsAsFactors = FALSE)
meta$X <- paste0("X", meta$X)
# Cell/sample name fix :
meta$sample[meta$sample == "sgSafe14_MT2179_2"] <- "sgSafe14_MT2179"

message("Reading counts matrix ...")
counts <- read.csv(file.path(data_dir, "kmt2a_filtered_counts_matrix.csv"), header = TRUE, stringsAsFactors = FALSE)
gene_ids <- counts$X
counts <- counts[, -1]
rownames(counts) <- gene_ids

# --- Ensure cell ID alignment & duplicates ---
if (!identical(sort(colnames(counts)), sort(meta$X))) {
  stop("Mismatch between cell IDs in counts and metadata! Please check.")
}
stopifnot(!anyDuplicated(meta$X))
stopifnot(!anyDuplicated(colnames(counts)))
meta <- meta[match(colnames(counts), meta$X), ]
stopifnot(identical(meta$X, colnames(counts)))
message("Cell ID alignment confirmed.")

# --- SCE object ---
message("Creating SingleCellExperiment object...")
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)))
colData(sce) <- DataFrame(meta)

# --- Scanpy-compatible CPM normalization ---
message("Normalizing with CPM + log1p (natural log, Scanpy compatible)...")
lib_sizes <- colSums(counts(sce))
lib_size_summary <- lib_sizes  # for output
CPM <- sweep(counts(sce), 2, lib_sizes / 1e6, "/")
logcounts(sce) <- log(CPM + 1)  # natural log, matches Scanpy sc.pp.log1p
message("Normalization (Scanpy-style) done.")

rowData(sce)$primerid <- rownames(sce)
colData(sce)$wellKey <- colnames(sce)
colData(sce)$cngeneson <- scale(colSums(counts(sce) > 0)) # cngeneson covariate

# --- Factor levels ---
expected_groups <- c("unenriched", "sgKmt2a_enriched")
colData(sce)$Groups <- factor(colData(sce)$Groups, levels = expected_groups)
colData(sce)$Groups <- relevel(colData(sce)$Groups, ref = "unenriched")
colData(sce)$sample <- factor(colData(sce)$sample)

message("Cells per group:")
print(table(colData(sce)$Groups))
message("Cells per sample:")
print(table(colData(sce)$sample))

# --- Stricter gene filtering: per sample/cell ---
MIN_CELLS_PER_SAMPLE <- 5
MIN_SAMPLES <- 3

message("Filtering genes: must be expressed in >=", MIN_CELLS_PER_SAMPLE, " cells in >=", MIN_SAMPLES, " samples ...")
samples <- levels(colData(sce)$sample)
expr_by_sample <- sapply(samples, function(s) {
  rowSums(counts(sce)[, colData(sce)$sample == s] > 0)
})
genes_pass <- rowSums(expr_by_sample >= MIN_CELLS_PER_SAMPLE) >= MIN_SAMPLES
genes_ok <- names(which(genes_pass))

message("Keeping ", length(genes_ok))
sce <- sce[genes_ok, ]

# --- Gene annotation filtering (exclude unwanted types) ---
exclude_types <- c(
  "pseudogene", "gene segment", "unclassified gene", "polymorphic pseudogene",
  "QTL", "transgene", "DNA segment", "unclassified other genome feature",
  "CTCF binding site", "intronic regulatory region", "imprinting control region",
  "enhancer", "complex/cluster/region", "transcription factor binding site"
)
message("Reading gene annotation ...")
gene_info_file <- file.path(data_dir, "MGIBatchReport_20250815_192205.txt")
if (!file.exists(gene_info_file)) stop("Gene annotation file not found: ", gene_info_file)
gene_info <- read.delim(gene_info_file, stringsAsFactors = FALSE)
message("The number of genes in gene_info_file is", nrow(gene_info))
gene_info_incl <- gene_info %>% filter(!Feature.Type %in% exclude_types)

genes_good <- intersect(rownames(sce), gene_info_incl$Input)
sce <- sce[genes_good, ]
message("After annotation filter: keeping ", length(genes_good))

# --- Build SingleCellAssay (recreate with filtered genes) ---
message("Recreating SingleCellAssay with filtered genes...")
sca <- FromMatrix(exprsArray = logcounts(sce),
                  cData   = colData(sce),
                  fData   = rowData(sce))

# --- Parallel MAST on blocks of genes ---
N_CORES <- 36
bp <- MulticoreParam(workers = N_CORES, RNGseed = 12345)

genes_all <- rownames(sce)
gene_blocks <- split(genes_all, cut(seq_along(genes_all), N_CORES, labels = FALSE))

mast_on_block <- function(this_genes) {
  message("Processing block with ", length(this_genes), " genes. First few: ",
          paste(head(this_genes, 5), collapse = ", "), if (length(this_genes) > 5) "..." else "")
  tryCatch({
    sca_block <- sca[this_genes, ]
    zlm_block <- zlm(
      formula = ~ Groups + cngeneson + (1 | sample),
      sca = sca_block,
      method = "glmer",
      ebayes = FALSE,
      strictConvergence = FALSE,
      fitArgsD = list(nAGQ = 0)
    )
    contrast_name <- "GroupssgKmt2a_enriched"
    summary_block <- summary(zlm_block, doLRT = contrast_name)
    summaryDt_block <- summary_block$datatable
    fcHurdle_block <- merge(
      summaryDt_block[contrast == contrast_name & component == "H", .(primerid, `Pr(>Chisq)`)],
      summaryDt_block[contrast == contrast_name & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
      by = "primerid"
    )
    # Only report log2 fold change
    fcHurdle_block[, log2FC := coef / log(2)] # Log2 fold change only
    fcHurdle_block$model_type <- "MAST_RE_with_cngeneson"
    return(fcHurdle_block)
  }, error = function(e) {
    message("Error in block with genes: ",
            paste(head(this_genes, 5), collapse=", "), if(length(this_genes)>5) "..." else "",
            "\nError: ", conditionMessage(e))
    return(data.table())
  })
}

message("Running parallel hurdle model (MAST) on ", N_CORES, " gene blocks ...")
results_list <- bplapply(gene_blocks, mast_on_block, BPPARAM=bp)
results_list <- Filter(function(d) !is.null(d) && nrow(d) > 0, results_list) # drop empty tables if any

fcHurdle_full <- rbindlist(results_list)
fcHurdle_full[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]  # global FDR
sig_genes_full <- fcHurdle_full[fdr < 0.05, , drop = FALSE]

message("Significant genes WITH cngeneson: ", nrow(sig_genes_full))

# --- Output CSVs ---
write.csv(fcHurdle_full, 
          file = file.path(output_dir, "sgKmt2a_MAST_DE_results.csv"), row.names = FALSE)

# --- Save all R objects ---
save(sca,
     fcHurdle_full, sig_genes_full,
     lib_size_summary,
     file = file.path(output_dir, "Kmt2a_MAST_analysis_objects_para.RData"),
     compress=TRUE, version=2
)

# --- Save session info ---
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sgKmt2a_session_info_para.txt"))

message("Analysis finished successfully for model WITH cngeneson (parallelized, block-wise).")