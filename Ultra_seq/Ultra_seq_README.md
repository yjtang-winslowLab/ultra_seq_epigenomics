# Ultra-Seq Analysis Pipeline

**Ultra-Seq: Massively Parallel Quantification of Growth Phenotype by Sequencing**

Created by **Yuning J. Tang**  
Dept of Genetics, Stanford University

---

## Overview

Tuba-seqUltra (a.k.a Ultra-Seq) analyzes growth phenotypes with gene perturbation from barcode sequencing data. The pipeline processes sample count files, calculates cell numbers with spike-in normalization, bootstraps tumor metrics, and performs statistical analysis to produce robust gene- and guide-level results.

---

## Pipeline Structure
Step 1: Data Merge      ──> Step 2: Spike-in Normalization ──> Step 3: Nested Dictionary  
       │                                      │
       └─────────────────┐                    │
                         │                    │
         Step 4: Tumor Bootstrap Metrics  ──>  Step 5: Serial Bootstrap      
                                              │
                                         Step 6: Stats, CIs, p-values, FDR, Meta Analysis

---

## Usage

Run each step sequentially. For serial bootstrapping (a.k.a. two-step bootstrap), **it is highly recommended that you use the script with parallelization.** For familiarizing yourself with the pipeline or working with very small sample size, you can use the single-loop script.

A sample dataset is included in the pilot_data folder. 

### 1. Prepare Your Data

* Place your sample count files in the appropriate data directory.
* Prepare an experimental information file (e.g. `expt_info.csv`).

### 2. Run Each Step

1. **Step 1:** Merge Data  
2. **Step 2:** Spike-in Normalization  
3. **Step 3:** Create Tumor Cell Number Dictionary  
4. **Step 4:** Bootstrap Tumor Metrics  
5. **Step 5:** Serial Bootstrap (metrics over many iterations)  
6. **Step 6:** Statistical Summary, CIs, FDR, Gene-Level Meta Analysis

---

## Citation

If you use this pipeline, please cite: Tang YJ et al.
