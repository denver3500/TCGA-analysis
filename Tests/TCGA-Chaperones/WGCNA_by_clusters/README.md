# WGCNA by Patient Clusters Analysis Pipeline

This directory contains scripts for performing WGCNA (Weighted Gene Co-expression Network Analysis) on TCGA data, stratified by patient clusters identified from previous chaperone expression analysis.

## Directory Structure
```
WGCNA_by_clusters/
├── 01_Normalization_by_clusters.R       # Data normalization for each patient cluster
├── 02_Network_construction.R            # WGCNA network construction
├── 03_Chaperone_module_analysis.R       # Analyze chaperone gene modules
├── 04_Cross_cluster_comparison.R        # Compare modules across clusters
├── processed_data/                      # Normalized data for each cluster
├── results/                             # WGCNA results and networks
├── chaperone_analysis/                  # Chaperone-specific analysis
└── cross_cluster_comparison/            # Cross-cluster comparison results
```

## Pipeline Overview

### 1. Data Normalization (`01_Normalization_by_clusters.R`)
- **Input**: Patient cluster assignments from `deseq2_clustering/`
- **Process**: 
  - Loads patient cluster data (e.g., TCGA-BLCA has clusters 1 and 2)
  - For each cluster, extracts genome-wide expression data
  - Applies DESeq2 normalization and variance stabilizing transformation
  - Prepares data in WGCNA format (samples as rows, genes as columns)
- **Output**: Normalized expression matrices for each cluster

### 2. Network Construction (`02_Network_construction.R`)
- **Input**: Normalized expression data from step 1
- **Process**:
  - Determines optimal soft-thresholding power for each cluster
  - Constructs co-expression networks using WGCNA
  - Identifies gene modules using hierarchical clustering
  - Calculates module eigengenes
- **Output**: 
  - Gene-to-module assignments
  - Module eigengenes
  - Network topology plots
  - Dendrogram visualizations

### 3. Chaperone Module Analysis (`03_Chaperone_module_analysis.R`)
- **Input**: Gene modules from step 2, chaperone gene list
- **Process**:
  - Maps chaperone genes to identified modules
  - Analyzes distribution of chaperone genes across modules
  - Identifies dominant modules containing chaperone genes
- **Output**:
  - Chaperone gene module assignments
  - Distribution plots and statistics
  - Module enrichment analysis

### 4. Cross-Cluster Comparison (`04_Cross_cluster_comparison.R`)
- **Input**: Module assignments from all clusters
- **Process**:
  - Compares chaperone module assignments between clusters within projects
  - Compares same cluster types across different cancer types
  - Calculates consistency metrics
- **Output**:
  - Comparison matrices
  - Consistency statistics
  - Cross-project heatmaps

## Key Features

### Cluster-Based Analysis
- **Patient Stratification**: Uses pre-computed patient clusters instead of expression-based stratification
- **Separate Networks**: Builds independent co-expression networks for each patient cluster
- **Comparative Analysis**: Identifies cluster-specific vs. shared molecular modules

### Quality Control
- **Sample Filtering**: Requires minimum 10 samples per cluster
- **Gene Filtering**: Removes low-expressed genes (≥10 counts in ≥25% samples)
- **Outlier Detection**: Identifies and handles outlier samples

### Biological Focus
- **Chaperone Genes**: Specifically tracks chaperone gene behavior across clusters
- **Module Consistency**: Measures how consistently chaperone genes cluster together
- **Cross-Cancer Comparison**: Identifies universal vs. cancer-specific patterns

## Usage

Run scripts in order:

```bash
# 1. Normalize data by clusters
Rscript 01_Normalization_by_clusters.R

# 2. Construct co-expression networks
Rscript 02_Network_construction.R

# 3. Analyze chaperone gene modules
Rscript 03_Chaperone_module_analysis.R

# 4. Compare across clusters
Rscript 04_Cross_cluster_comparison.R
```

## Expected Results

### For Each Cancer Type:
- **Cluster 1 Network**: Co-expression modules for patient cluster 1
- **Cluster 2 Network**: Co-expression modules for patient cluster 2
- **Comparison**: How chaperone genes behave differently between clusters

### Cross-Cancer Analysis:
- **Cluster 1 Comparison**: How cluster 1 networks compare across cancer types
- **Cluster 2 Comparison**: How cluster 2 networks compare across cancer types
- **Universal Patterns**: Chaperone modules consistent across cancers

## Biological Interpretation

### Cluster-Specific Modules:
- **Unique Networks**: Each patient cluster may have distinct co-expression patterns
- **Therapeutic Implications**: Different clusters may respond to different treatments
- **Biomarker Discovery**: Cluster-specific modules could serve as biomarkers

### Chaperone Gene Behavior:
- **Co-expression Patterns**: Which chaperone genes work together in each cluster
- **Functional Modules**: Groups of chaperone genes with related functions
- **Cancer-Specific Roles**: How chaperone networks differ between cancer types

## Technical Details

### WGCNA Parameters:
- **Soft Power**: Automatically selected based on scale-free topology
- **Module Size**: Minimum 30 genes per module
- **Merge Threshold**: 0.25 for merging similar modules
- **Deep Split**: 2 for sensitive module detection

### Statistical Methods:
- **DESeq2**: For proper RNA-seq normalization
- **VST**: Variance stabilizing transformation
- **Hierarchical Clustering**: Ward's method for module detection
- **Correlation**: Pearson correlation for co-expression

This pipeline provides comprehensive analysis of how patient molecular subtypes influence gene co-expression networks, with specific focus on chaperone gene behavior.
