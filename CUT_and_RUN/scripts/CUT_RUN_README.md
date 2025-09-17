# CUT&RUN Pipeline

## Steps

1. **Merge FASTQ batches**: Combine raw FASTQ.gz files per sample.
2. **QC with FastQC**: Automated quality control for merged FASTQ files.
3. **Trim adapters with Trimmomatic**: Batch trims paired-end reads.
4. **Bowtie2 alignment (genome)**: Align trimmed reads to reference genome.
5. **Spike-in alignment/control**: Align to *E. coli* for spike-in normalization; compute sequencing depth.
6. **SAM→Fragment BED**: Converts alignments to fragment-format BEDs (for coverage & counting).
7. **Fragment binning/counts**: Counts fragments per genomic bin (e.g., 500 bp).
8. **Normalized bedGraph**: Generates scaled coverage tracks per sample.
9. **Peak calling with SEACR**: Calls peaks vs control for each target protein.
10. **BAM→bigWig**: Sort/index BAM, produce bigWig tracks for visualization.
11. **deepTools matrix/plots**: Generates heatmaps/profiles of signal over regions.
12. **Consensus peaks**: Intersects replicate peaks for consensus BED.
13. **Peak annotation (HOMER)**: Annotate peaks with genomic features.
14. **Scaled/smoothed bigWigs**: Final bigWigs for publication, normalized and smoothed.

---

## Folder Structure

project/
├── raw_data/               # Contains CR_1 ... CR_N folders of raw FASTQ.gz
├── fastqc_output/
├── trimmed_data/
├── alignment/
│   └── sam/
├── analysis/
│   ├── bed/
│   ├── bedgraph/
│   ├── matrix/
│   └── plots/
├── peakCalling/
│   └── SEACR/
├── bam/
├── bigwig/
├── bigwigs_scaled_smoothed/
└── logs/

## Requirements
FastQC
Trimmomatic
Bowtie2
Samtools
Bedtools
deepTools
SEACR
HOMER

## Reference files
Bowtie2 index for your genome (e.g., mm10)
Adapter sequences file for Trimmomatic
Chromosome sizes file (for bedtools genomecov)
Blacklist file (optional, for deepTools)
Region annotation file (BED, for deepTools plots)

## Usage
General usage:
Run each script in sequence, editing the paths and variables at the top to match your directories, sample names, and parameters.
Note: scaling factors for normalization need to be calculated in R using QC_NormScale.R

Example:
bash CR_step01_merge_fastq_batches.sh
bash CR_step02_batch_fastqc.sh
bash CR_step03_batch_trimmomatic.sh
...
bash CR_step14_make_scaled_smoothed_bigwig.sh
Scripts can be run manually, or included in your own workflow manager (snakemake/nextflow). Each script produces logs and informative output for troubleshooting.

## Credits
Written by Y.J. Tang, 

Original scripts: Designed for mouse lung adenocarcinoma cell line CUT&RUN data.

License
MIT License