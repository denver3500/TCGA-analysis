# Heatmap Publication Figure Guide

## Available Heatmap Types

### 1. Correlation Heatmaps
- **Location**: `TCGA-Chaperones/pictures/correlations_pheatmap/`
- **Files**: `heatmap_TCGA-*.png/pdf`
- **Description**: Gene-gene correlation matrices for chaperone genes in each cancer type

### 2. Patient Stratification Heatmaps
- **Location**: `TCGA-Chaperones/heatmap_patients_stratification/deseq2_clustering/*/`
- **Files**: `*_bidirectional_clustered_heatmap.pdf`
- **Description**: Patient clustering based on chaperone expression patterns

### 3. Pathway Enrichment Heatmaps
- **Location**: `TCGA-Chaperones/heatmap_patients_stratification/`
- **Files**: `pathways_heatmap.pdf`, `top_enriched_pathways.pdf`
- **Description**: Functional pathway analysis results

## Generated Publication Figures

### Figure 1: Multi-cancer correlation heatmaps
- **File**: `Figure_1_Multi_cancer_correlation_heatmaps.png/pdf`
- **Description**: 4x2 grid showing chaperone gene correlation patterns across 8 cancer types
- **Use for**: Demonstrating consistency/variability of chaperone co-expression patterns

### Figure 2: Patient stratification heatmaps
- **File**: `Figure_2_Patient_stratification_heatmaps.png/pdf`
- **Description**: 3x2 grid showing patient clustering based on chaperone expression
- **Use for**: Showing how chaperone expression can stratify patients

### Figure 3: Comprehensive heatmap workflow
- **File**: `Figure_3_Comprehensive_heatmap_workflow.png/pdf`
- **Description**: Multi-panel figure showing different types of heatmap analyses
- **Use for**: Main figure demonstrating the heatmap-based analytical approach

### Figure 4: Cancer type heatmap comparison
- **File**: `Figure_4_Cancer_type_comparison.png/pdf`
- **Description**: Side-by-side comparison of correlation heatmaps
- **Use for**: Highlighting specific differences between cancer types

## Customization Options

### Available Cancer Types for Correlation Heatmaps:
- TCGA-BRCA (Breast Cancer)
- TCGA-CHOL (Cholangiocarcinoma)
- TCGA-COAD (Colon Adenocarcinoma)
- TCGA-ESCA (Esophageal Carcinoma)
- TCGA-HNSC (Head and Neck Squamous Cell Carcinoma)
- TCGA-LIHC (Liver Hepatocellular Carcinoma)
- TCGA-LUAD (Lung Adenocarcinoma)
- TCGA-LUSC (Lung Squamous Cell Carcinoma)
- TCGA-PRAD (Prostate Adenocarcinoma)
- TCGA-READ (Rectum Adenocarcinoma)
- TCGA-STAD (Stomach Adenocarcinoma)

### Available Cancer Types for Patient Stratification:
- TCGA-BLCA, TCGA-BRCA, TCGA-CESC, TCGA-CHOL, TCGA-COAD
- TCGA-ESCA, TCGA-HNSC, TCGA-KIRC, TCGA-KIRP, TCGA-LIHC
- TCGA-LUAD, TCGA-LUSC, TCGA-PAAD, TCGA-PRAD, TCGA-READ
- TCGA-SARC, TCGA-STAD, TCGA-THYM

## Usage Tips

1. **Focus on specific cancer types**: Modify the `selected_cancers` vector to include only the cancer types relevant to your study
2. **Adjust grid layouts**: Change `ncol` and `nrow` parameters in `wrap_plots()` to customize layouts
3. **Modify titles and labels**: Update `plot_annotation()` titles and subtitles
4. **Change figure sizes**: Adjust `width` and `height` parameters in `ggsave()`
5. **Add annotations**: Use `plot_annotation()` to add figure labels (A, B, C, etc.)

## File Format Information

- **PDF files**: Vector graphics, best for publications
- **PNG files**: Raster graphics, good for presentations
- **High resolution**: All figures saved at 300 DPI for publication quality

