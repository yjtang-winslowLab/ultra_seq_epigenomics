# RNA-Seq Analysis Pipeline

This repository contains shell scripts for bulk RNA-seq analysis of Kat7-deficient vs. Kat7-proficient AT2 cells.  
The workflow covers quality control, alignment, and quantification using FastQC, STAR, and Subread (featureCounts).

> **Note:** The scripts assume that the mm39 genome has already been indexed for STAR.  
> For genome indexing instructions, see the WM3835 scripts included in this repository.

---

## Script Overview

- **run_fastqc.sh**: Runs FastQC for quality control on raw FASTQ files.
- **star_align.sh**: Aligns paired-end RNA-seq samples to the reference genome using STAR.
- **featurecounts.sh**: Performs read quantification using featureCounts and a GTF annotation.

---

## Requirements

- A UNIX-like environment (Linux/macOS)
- [STAR](https://github.com/alexdobin/STAR)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Subread (featureCounts)](http://subread.sourceforge.net/)
- Ability to load tools via environment modules (or update the scripts for your system)

---

## Setup Instructions

**Before running the scripts, edit path variables at the top of each script to match your system:**
- Replace example paths (`/your/reference/...`) with the locations of your reference files and sample data.
- Update sample names or sample lists if needed.

---

## Usage

### 1. FASTQ Quality Control

Run:
```bash
bash run_fastqc.sh
```
Results are written to ${base_dir}/fastqc/.

### 2. STAR Alignment
Edit base_dir and index in star_align.sh.
Run:
```bash
bash star_align.sh
```
Aligned BAM files are written to ${base_dir}/alignment/star_align/.

### 3. featureCounts Quantification
Run:
```bash
bash featurecounts.sh
```
Produces a count matrix (featureCounts_results.txt) in ${base_dir}/counts/.

### 4. General Notes
Scripts expect paired-end FASTQ files named as SAMPLE_1.fq.gz and SAMPLE_2.fq.gz in your data directory.
If you use a different environment setup (e.g., not using modules), adjust or remove the module load lines.
Output directories will be created as needed.
Make scripts executable (chmod +x scriptname.sh) if desired.

License
Distributed under the MIT License.