# RNA-Seq Analysis Pipeline

This repository contains scripts for bulk RNA-seq analysis workflows for KP lung adenocarcinoma cell lines treated with 1uM WM3835 or DMSO in vitro. 

The shell scripts include quality control, alignment, and quantification. Tools included: FastQC, STAR, and Subread (featureCounts).

## Scripts Overview

- **star_genome_generate.sh**: Generates the STAR genome index using a reference FASTA and GTF.
- **run_fastqc.sh**: Runs FastQC for quality control on raw FASTQ files.
- **star_align.sh**: Aligns paired-end RNA-seq samples to the reference genome with STAR.
- **featurecounts.sh**: Counts reads mapping to genes using featureCounts and a GTF annotation.

---

## Requirements

- A UNIX-like environment (Linux/macOS)
- [STAR](https://github.com/alexdobin/STAR)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Subread (featureCounts)](http://subread.sourceforge.net/)
- Environment modules (or load tools by your preferred method)

---

## Setup Instructions

**Edit each script to set your own file system paths:**
- Replace `/your/reference/...` or other example folders with actual locations for reference files and data on your system.
- Update lists of sample names or BAM files as appropriate.

---

## Usage

### 1. STAR Genome Indexing

Edit by setting the reference FASTA and GTF paths at the top of `star_genome_generate.sh`.
Run:
```bash
bash star_genome_generate.sh
```

### 2. FASTQ Quality Control
Run:
```bash
bash run_fastqc.sh
```
Results are written to ${base_dir}/fastqc/.

### 3. STAR Alignment
Edit base_dir and index in star_align.sh.
Run:
```bash
bash star_align.sh
```
Aligned BAM files are written to ${base_dir}/alignment/star_align/.

### 4. featureCounts Quantification
Run:
```bash
bash featurecounts.sh
```
Produces a count matrix (featureCounts_results.txt) in ${base_dir}/counts/.

### 5. General Notes
Scripts expect paired-end FASTQ files named as SAMPLE_1.fq.gz and SAMPLE_2.fq.gz in your data directory.
If you use a different environment setup (e.g., not using modules), adjust or remove the module load lines.
Output directories will be created as needed.
Make scripts executable (chmod +x scriptname.sh) if desired.

License
Distributed under the MIT License.