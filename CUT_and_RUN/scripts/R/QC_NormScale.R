# CUT&RUN Analysis

# --- Setup ---
library(tidyverse)
library(viridis)
library(ggpubr)
library(corrplot)

# --- User Parameters ---
sampleList <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8")
targetList <- c("IgG", "H3K4me3", "H3K14ac", "H3K14ac", "MLL1", "MLL1", "KAT7", "KAT7")

# DIR PATHS (edit if needed, but best to keep relative for GitHub)
bowtie2_summary_path <- "data/bowtie2_summary/"
spike_summary_path <- "data/bowtie2_Spikein_summary/"
bed_path <- "data/bin_bed/"
outdir <- "results/"
if(!dir.exists(outdir)) dir.create(outdir)

# --- 1. Mapping Summary (main genome) ---
alignResults <- list()
for (sample in sampleList) {
  fn <- file.path(bowtie2_summary_path, paste0(sample, "_bowtie2.txt"))
  if(!file.exists(fn)) next
  lines <- readLines(fn)
  lines <- lines[!grepl("^Warning:", lines)]
  df <- read.table(text = lines, header = FALSE, fill = TRUE)
  alignRate <- substr(as.character(df$V1[6]), 1, nchar(as.character(df$V1[6])) - 1)
  alignResults[[sample]] <- data.frame(
    Sample_id = sample,
    SequencingDepth = as.numeric(df$V1[1]),
    MappedFragNum_mm10 = as.numeric(df$V1[4]) + as.numeric(df$V1[5]),
    AlignmentRate_mm10 = as.numeric(alignRate)
  )
}
alignResult <- bind_rows(alignResults)
alignResult <- alignResult %>% mutate(AlignmentRate_mm10 = paste0(AlignmentRate_mm10, "%"))

# --- 2. Mapping Summary (spike-in) ---
spikeResults <- list()
for (sample in sampleList) {
  fn <- file.path(spike_summary_path, paste0(sample, "_bowtie2_spikeIn.txt"))
  if(!file.exists(fn)) next
  lines <- readLines(fn)
  lines <- lines[!grepl("^Warning:", lines)]
  df <- read.table(text = lines, header = FALSE, fill = TRUE)
  alignRate <- substr(as.character(df$V1[6]), 1, nchar(as.character(df$V1[6])) - 1)
  spikeResults[[sample]] <- data.frame(
    Sample_id = sample,
    SequencingDepth = as.numeric(df$V1[1]),
    MappedFragNum_spikeIn = as.numeric(df$V1[4]) + as.numeric(df$V1[5]),
    AlignmentRate_spikeIn = as.numeric(alignRate)
  )
}
alignResult_spikeIn <- bind_rows(spikeResults)
alignResult_spikeIn <- alignResult_spikeIn %>%
  mutate(AlignmentRate_spikeIn_pct = paste0(AlignmentRate_spikeIn, "%"))

# --- 3. Summarize Both Alignments ---
alignSummary <- left_join(alignResult, alignResult_spikeIn, by = c("Sample_id", "SequencingDepth"))
alignSummary$Target <- targetList

write.csv(alignSummary, file.path(outdir, "alignSummary.csv"), row.names = FALSE)
print(alignSummary)

# --- 4. QC Plots ---
p1 <- ggplot(alignSummary, aes(Target, SequencingDepth/1e6, fill=Target)) +
  geom_boxplot() +
  geom_jitter(aes(color=Target), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete=TRUE, option="magma") +
  scale_color_viridis(discrete=TRUE) +
  theme_bw() + ylab("Sequencing Depth (Million)") + ggtitle("A. Sequencing Depth")

p2 <- ggplot(alignSummary, aes(Target, MappedFragNum_mm10/1e6, fill=Target)) +
  geom_boxplot() +
  geom_jitter(aes(color=Target), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete=TRUE, option="magma") +
  scale_color_viridis(discrete=TRUE) +
  theme_bw() + ylab("Mapped Fragments (Million, mm10)") + ggtitle("B. Alignable Fragments")

p3 <- ggplot(alignSummary, aes(Target, as.numeric(gsub("%","",AlignmentRate_mm10)), fill=Target)) +
  geom_jitter(aes(color=Target), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete=TRUE, option="magma") +
  scale_color_viridis(discrete=TRUE) +
  theme_bw() + ylab("% of Mapped Frags (mm10)") + ggtitle("C. Alignment Rate (mm10)")

p4 <- ggplot(alignSummary, aes(Target, AlignmentRate_spikeIn, fill=Target)) +
  geom_jitter(aes(color=Target), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete=TRUE, option="magma") +
  scale_color_viridis(discrete=TRUE) +
  theme_bw() + ylab("Spike-in Alignment Rate") + ggtitle("D. Alignment Rate (spike-in)")

ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
ggsave(file.path(outdir, "all_QC_plots.png"), dpi=150, width=12, height=8)

# --- 5. Scaling Factor Calculation ---
scaleFactorOut <- alignResult_spikeIn %>%
  select(Sample_id, AlignmentRate_spikeIn) %>%
  mutate(scaleFactor = 1/(AlignmentRate_spikeIn))
write.csv(scaleFactorOut, file.path(outdir, "scaleFactors.csv"), row.names = FALSE)
print(scaleFactorOut)

alignSummary <- left_join(alignSummary, scaleFactorOut, by = "Sample_id")

# --- 6. Reproducibility Across Samples ---
fragCount <- NULL
for (sample in sampleList) {
  fn <- file.path(bed_path, paste0(sample, "_bowtie2.fragmentsCount.bin500.bed"))
  if(!file.exists(fn)) next
  fragCountTmp <- read.table(fn, header = FALSE)
  colnames(fragCountTmp) <- c("chrom", "bin", sample)
  if (is.null(fragCount)) {
    fragCount <- fragCountTmp
  } else {
    fragCount <- full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}

mat <- fragCount %>% select(-chrom, -bin)
mat_log2 <- log2(mat + 1) # Added +1 for log2 safety

colnames(mat_log2) <- sampleList
correlation_matrix <- cor(mat_log2, use = "complete.obs")
print(correlation_matrix)

png(file.path(outdir, "sample_reproducibility_correlation.png"), width=900, height=900)
corrplot(correlation_matrix, method="color", outline=TRUE, addgrid.col="darkgray", order="hclust",
         addrect=3, rect.col="black", rect.lwd=3, cl.pos="b",
         tl.col="indianred4", tl.cex=1, cl.cex=1,
         addCoef.col="black", number.digits=2, number.cex=1,
         col=colorRampPalette(c("midnightblue", "white", "darkred"))(100))
dev.off()