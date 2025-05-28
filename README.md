# scRNA-Seq Analysis Pipeline with Seurat

This repository contains an R script for the analysis of single-cell RNA sequencing (scRNA-Seq) data using the [Seurat](https://satijalab.org/seurat/) package. The dataset used includes peripheral blood mononuclear cells (PBMCs) from human samples (SRA: [SRA713577](https://www.ncbi.nlm.nih.gov/sra/?term=SRA713577)).

## Contents

- `scRNA-Seq.R`: The main script performing a full scRNA-Seq analysis pipeline including:
  - Data loading and preprocessing
  - Quality control
  - Normalization and scaling
  - Dimensionality reduction (PCA, UMAP, t-SNE)
  - Clustering and cell type annotation
  - Marker gene identification

## Requirements

Install the necessary R packages before running the script:

```r
install.packages("tidyverse")
install.packages("Seurat")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("dplyr")
# UMAP requires Python package
library(reticulate)
reticulate::py_install(packages = 'umap-learn')
```

## Usage
1. Prepare the dataset: download the following file and place it in the correct path:
```r
~/Downloads/trascrittomica/scRNA-seq/SRA713577_SRS3363004.sparse.RData
```
Make sure to update the load() line in the script if the file is located elsewhere.

2. Run the script: you can execute the entire pipeline by running the script in R or RStudio:
```r
source("scRNA-Seq.R")
```
The script will:
- Create a Seurat object from the sparse matrix
- Filter low-quality cells
- Normalize and scale the data
- Identify highly variable genes
- Perform PCA, UMAP, and t-SNE for dimensionality reduction
- Cluster the cells and identify marker genes
- Annotate cell types using canonical markers

## Cell Type Annotations
The final clusters are annotated with likely cell types:
Cluster	Assigned Cell Type
- Oligodendrocytes
- T-cells
- Granulocytes
- NK-cells
- B-cells
- Monocytes
- Excitatory neurons
- Langerhans cells
- Dendritic cells
- Adipocytes
- Basal respiratory cells
