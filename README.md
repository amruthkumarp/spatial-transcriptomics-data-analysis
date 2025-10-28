# 🧬 Spatial Transcriptomics Analysis of Human Prostate Cancer

## 📘 Overview
This repository contains code and data processing scripts for the **Spatial Transcriptomics (ST)** analysis of **Human Prostate Cancer** samples.  
The study explores spatial gene expression patterns across prostate tissue sections using **ST data from CNGBdb**.

---

## 📂 Repository Structure
data
---

## 🔗 Data Source
The raw spatial transcriptomics data used in this project was obtained from the **CNGBdb (China National GeneBank Database)**:

**Dataset Link:** [STDS0000131 – Human Prostate Cancer Spatial Transcriptomics](https://db.cngb.org/stomics/datasets/STDS0000131/summary)

> The dataset provides Visium spatial transcriptomics data from human prostate cancer tissue samples, including gene expression matrices and histological images.

---

## ⚙️ Analysis Workflow

### 1️⃣ Quality Control (`Qc.r`)
- Reads the raw `.h5ad` / matrix data.  
- Filters out low-quality spots based on:
  - Low total counts
  - Low gene numbers per spot
  - High mitochondrial gene percentages  
- Produces violin plots and scatter plots for QC metrics.  

### 2️⃣ Spatial Object Creation (`spatial_obj.r`)
- Constructs a **SpatialExperiment** or **Seurat** object.
- Integrates spatial coordinates with gene expression.
- Links tissue images to molecular data.
- Performs normalization and visualization setup for downstream spatial analysis.

### 3️⃣ Downstream Analysis (to be added)
Planned future analyses include:
- Spatial clustering
- Identification of spatially variable genes
- Differential expression across tumor and normal regions
- Visualization overlays (gene expression heatmaps on tissue)

---

## 🧰 Requirements

### R Dependencies
Install the following R packages before running the scripts:
```r
install.packages(c("Seurat", "SpatialExperiment", "ggplot2", "dplyr", "patchwork"))
