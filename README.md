# Single-cell-variability-project
# Early Transcriptomic Heterogeneity and Evolutionary Constraints in *Danio rerio*

## About The Project
This repository contains the computational pipeline developed to investigate the origins of gene expression variability and evolutionary constraints during early embryogenesis. Utilizing single-cell RNA sequencing (scRNA-seq) data from the *Danio rerio* (zebrafish) model at the 10 days post-fertilization (10 dpf) larval stage, this project explores whether the regulatory architecture governing tissue-specific transcriptional plasticity is established prior to macroscopic organ maturation.

This work was conducted at the **Department of Ecology and Evolution (DEE)** at the **University of Lausanne (UNIL)**.

## Scientific Objectives
*   **Isolate Biological Stochasticity:** Separate true biological variance from technical noise using advanced mean-variance modeling.
*   **Evaluate Transcriptomic Compartmentalization:** Map specific functional gene networks (e.g., developmental growth, tissue regeneration, sensory perception) across the spatial topology of the larval organism.
*   **Assess Pleiotropic Constraints:** Compare the biological variance of focal Highly Variable Genes (HVGs) from the Central Nervous System (CNS) against non-target peripheral tissues (liver, muscle, fin) to demonstrate transcriptomic canalization.

## Dataset
The analyses utilize raw single-cell transcriptomic data extracted from the public **Zebra Hub** database (Multimodal Zebrafish Developmental Atlas). 
*   **Format:** The raw data is initially provided in Python's AnnData format (`.h5ad`) and bridged into the R/Seurat ecosystem.
*   **Scale:** The global matrix lists the expression levels of 32,060 genes across 20,579 unique cells from 4 independent biological replicates.

## Pipeline Overview
The analytical workflow is structured into sequential modules:

1.  **Data Import & Initialization:** Conversion of `.h5ad` matrices to R structures via the `reticulate` package and initialization of the `Seurat` object.
2.  **Quality Control (QC) & Filtering:** Removal of empty droplets, doublets, and apoptotic/lysed cells based on unique feature counts (nFeature_RNA) and mitochondrial transcript ratios (pct_counts_mt).
3.  **Normalization & Variance Modeling:** Implementation of the `scran` algorithm utilizing cellular pre-pooling (deconvolution) to compute robust size factors, followed by local regression to isolate biological variance from technical background noise.
4.  **Dimensionality Reduction & Spatial Mapping:** Application of Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP) to visualize functional gene modules in a 2D latent space.
5.  **Functional Enrichment:** Gene Ontology (GO) analysis to identify significantly enriched Biological Processes among the top 10% isolated HVGs.
6.  **Inter-Tissue Canalization Analysis:** Targeted subsetting of specific tissues (CNS vs. periphery) coupled with independent variance modeling and non-parametric statistical testing (Wilcoxon rank-sum tests) to quantify the repression of gene expression variability outside of focal organs.

## Prerequisites and Dependencies
All bioinformatics analyses and data manipulations are carried out within the R programming environment (v4.5.3). Ensure the following packages are installed:

*   `Seurat` (v5.5.0) - Core single-cell analytical ecosystem
*   `scran` (v1.38.1) - Variance modeling and normalization
*   `SingleCellExperiment` - Data structural framework
*   `clusterProfiler` - Functional enrichment (Gene Ontology)
*   `org.Dr.eg.db` - Zebrafish genome annotation database
*   `reticulate` - Python/R interoperability
*   `tidyverse` / `ggplot2` / `ggpubr` - Data manipulation and statistical visualization

## Usage
Scripts should be executed in chronological order. Due to the high computational memory required for deconvolution and spatial mapping, it is recommended to run the pipeline on a machine with at least 16GB of RAM or a dedicated computing cluster. 

## Acknowledgments
*   **Author:** Nourine FAIZ (Assistant Engineer Intern, ESTBB - Catholic University of Lyon)
*   **Supervisor:** Smith SAM
*   **Host Institution:** University of Lausanne (UNIL), Department of Ecology and Evolution (DEE), Robinson-Rechavi Lab
# Bioinformatics pipeline summary diagram for scRNA-seq analysis.
<img width="991" height="619" alt="image" src="https://github.com/user-attachments/assets/e657173e-2d9d-4099-a424-4f6492261ccc" />
