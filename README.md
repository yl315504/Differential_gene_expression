# Differential Gene Expression Pipeline
## Overview
Differential gene expression (DGE) analysis is commonly performed as a downstream step in RNA-seq workflows. It can also be utilized in DNA-Protein-Crosslinking sequencing (DPC-seq).
This Bash script invokes an R script to carry out the analysis, employing **edgeR** as the primary tool.
The R script generates normalized count tables, diagnostic plots, and identifies differentially expressed genes for specified contrasts.

## Usage
```bash
Rscript edgeR.R <COUNTS_FILE> <SAMPLE_FILE> <OUT_DIR> <EDGER_BCV> <CONTRAST_FILE> <TPM_file> <AnnotationFile|none>

Parameters
COUNTS_FILE: Tab-delimited file with raw count matrix (genes x samples).
SAMPLE_FILE: CSV file mapping samples to experimental groups.
OUT_DIR: Directory to store output files.
EDGER_BCV: Biological coefficient of variation.
CONTRAST_FILE: CSV file specifying contrasts for group comparisons.
TPM_file: File containing transcript-per-million values.
AnnotationFile: Additional annotation file or "none".

Outputs
Normalized count tables (CPM, TPM).
Differential expression results for each contrast (logFC, p-value, FDR).
Diagnostic plots:
MA plots
Volcano plots
PCA plots
Heatmaps

Dependencies
R Packages: edgeR, DESeq2, gplots, RColorBrewer, NMF
Input Files: Properly formatted counts and sample files.
