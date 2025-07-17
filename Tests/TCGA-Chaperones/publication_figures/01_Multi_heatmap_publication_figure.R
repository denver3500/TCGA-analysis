#!/usr/bin/env Rscript

# Multi-Heatmap Publication Figure Generator
# This script creates publication-quality multi-panel figures combining various heatmaps

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(patchwork)
library(png)
library(jpeg)
library(magick)
library(pheatmap)
library(ComplexHeatmap)
library(readr)
library(dplyr)

# Create output directory
output_dir <- "TCGA-Chaperones/publication_figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
log_file <- file.path(output_dir, "multi_heatmap_log.txt")
con <- file(log_file, "w")
sink(con, split = TRUE)

message("=== Multi-Heatmap Publication Figure Generator ===")
message("Starting at: ", Sys.time())

# Define file paths for different types of heatmaps (focusing only on existing heatmaps)
heatmap_paths <- list(
  correlation_heatmaps = "TCGA-Chaperones/pictures/correlations_pheatmap",
  patient_stratification = "TCGA-Chaperones/heatmap_patients_stratification/deseq2_clustering",
  pathways = "TCGA-Chaperones/heatmap_patients_stratification"
)

# Function to read and prepare image for plotting
read_image_for_plot <- function(file_path, format = "auto") {
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(NULL)
  }
  
  tryCatch({
    if (format == "auto") {
      format <- tolower(tools::file_ext(file_path))
    }
    
    if (format %in% c("png", "PNG")) {
      img <- png::readPNG(file_path)
    } else if (format %in% c("jpg", "jpeg", "JPG", "JPEG")) {
      img <- jpeg::readJPEG(file_path)
    } else if (format %in% c("pdf", "PDF")) {
      # Convert PDF to PNG using magick
      img <- magick::image_read_pdf(file_path, density = 300)
      img <- magick::image_convert(img, "png")
      
      # Convert to array format that can be used with grid::rasterGrob
      img_array <- magick::image_data(img, "rgba")
      
      # Convert to proper format for R plotting
      img_dims <- dim(img_array)
      img_raster <- array(as.numeric(img_array), dim = img_dims)
      img_raster <- aperm(img_raster, c(2, 1, 3))  # transpose
      img_raster <- img_raster[, ncol(img_raster):1, ]  # flip vertically
      
      # Normalize to 0-1 range
      if (max(img_raster) > 1) {
        img_raster <- img_raster / 255
      }
      
      # Remove alpha channel if present
      if (dim(img_raster)[3] == 4) {
        img_raster <- img_raster[,,1:3]
      }
      
      return(img_raster)
    } else {
      # Try to convert using magick as fallback
      img <- magick::image_read(file_path)
      img <- magick::image_convert(img, "png")
      
      # Convert to array format
      img_array <- magick::image_data(img, "rgba")
      img_dims <- dim(img_array)
      img_raster <- array(as.numeric(img_array), dim = img_dims)
      img_raster <- aperm(img_raster, c(2, 1, 3))
      img_raster <- img_raster[, ncol(img_raster):1, ]
      
      if (max(img_raster) > 1) {
        img_raster <- img_raster / 255
      }
      
      if (dim(img_raster)[3] == 4) {
        img_raster <- img_raster[,,1:3]
      }
      
      return(img_raster)
    }
    
    return(img)
    
  }, error = function(e) {
    message("Error reading image: ", file_path, " - ", e$message)
    return(NULL)
  })
}

# Function to create a grid plot from image
create_image_plot <- function(img, title = "", subtitle = "") {
  if (is.null(img)) return(NULL)
  
  g <- grid::rasterGrob(img, interpolate = TRUE)
  
  p <- ggplot() + 
    annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() +
    labs(title = title, subtitle = subtitle) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

# =============================================================================
# APPROACH 1: Multi-cancer correlation heatmaps grid
# =============================================================================

message("\n--- Creating multi-cancer correlation heatmaps figure ---")

# Get all correlation heatmap files
correlation_files <- list.files(heatmap_paths$correlation_heatmaps, 
                               pattern = "heatmap_TCGA.*\\.png$", 
                               full.names = TRUE)

# Select representative cancer types for publication
selected_cancers <- c("TCGA-BRCA", "TCGA-LUAD", "TCGA-COAD", "TCGA-PRAD", 
                     "TCGA-HNSC", "TCGA-LIHC", "TCGA-STAD", "TCGA-READ")

correlation_plots <- list()
for (cancer in selected_cancers) {
  cancer_file <- correlation_files[grep(cancer, correlation_files)]
  if (length(cancer_file) > 0) {
    img <- read_image_for_plot(cancer_file[1])
    if (!is.null(img)) {
      plot_title <- gsub("TCGA-", "", cancer)
      correlation_plots[[cancer]] <- create_image_plot(img, title = plot_title)
    }
  }
}

if (length(correlation_plots) > 0) {
  # Create multi-panel figure
  multi_corr_plot <- wrap_plots(correlation_plots, ncol = 4, nrow = 2)
  
  # Add overall title
  multi_corr_plot <- multi_corr_plot + 
    plot_annotation(
      title = "Chaperone Gene Correlation Patterns Across Cancer Types",
      subtitle = "Correlation heatmaps showing co-expression patterns of chaperone genes",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )
  
  # Save the figure
  ggsave(file.path(output_dir, "Figure_1_Multi_cancer_correlation_heatmaps.png"), 
         multi_corr_plot, width = 16, height = 8, dpi = 300, bg = "white")
  
  ggsave(file.path(output_dir, "Figure_1_Multi_cancer_correlation_heatmaps.pdf"), 
         multi_corr_plot, width = 16, height = 8, bg = "white")
  
  message("Created multi-cancer correlation heatmaps figure")
}

# =============================================================================
# APPROACH 2: Patient stratification comparison
# =============================================================================

message("\n--- Creating patient stratification comparison figure ---")

# Get patient stratification heatmaps - prefer PNG files over PDF
stratification_files <- list()
cancer_types <- list.dirs(heatmap_paths$patient_stratification, full.names = FALSE, recursive = FALSE)
cancer_types <- cancer_types[grepl("^TCGA-", cancer_types)]

for (cancer in cancer_types) {
  # First try PNG file
  png_file <- file.path(heatmap_paths$patient_stratification, cancer, 
                       paste0(cancer, "_bidirectional_clustered_heatmap.png"))
  pdf_file <- file.path(heatmap_paths$patient_stratification, cancer, 
                       paste0(cancer, "_bidirectional_clustered_heatmap.pdf"))
  
  if (file.exists(png_file)) {
    stratification_files[[cancer]] <- png_file
  } else if (file.exists(pdf_file)) {
    stratification_files[[cancer]] <- pdf_file
  }
}

# Select first 6 for display
selected_stratification <- head(stratification_files, 6)
stratification_plots <- list()

for (cancer_type in names(selected_stratification)) {
  file_path <- selected_stratification[[cancer_type]]
  cancer_name <- gsub("TCGA-", "", cancer_type)
  
  img <- read_image_for_plot(file_path)
  if (!is.null(img)) {
    stratification_plots[[cancer_type]] <- create_image_plot(
      img, 
      title = paste(cancer_name, "Patient Clusters"),
      subtitle = "Bidirectional clustering"
    )
  }
}

if (length(stratification_plots) > 0) {
  multi_strat_plot <- wrap_plots(stratification_plots, ncol = 3, nrow = 2)
  
  multi_strat_plot <- multi_strat_plot + 
    plot_annotation(
      title = "Patient Stratification Based on Chaperone Expression",
      subtitle = "Bidirectional clustering reveals distinct patient subgroups",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )
  
  ggsave(file.path(output_dir, "Figure_2_Patient_stratification_heatmaps.png"), 
         multi_strat_plot, width = 18, height = 12, dpi = 300, bg = "white")
  
  ggsave(file.path(output_dir, "Figure_2_Patient_stratification_heatmaps.pdf"), 
         multi_strat_plot, width = 18, height = 12, bg = "white")
  
  message("Created patient stratification comparison figure")
}

# =============================================================================
# APPROACH 3: Comprehensive analysis workflow figure (heatmaps only)
# =============================================================================

message("\n--- Creating comprehensive heatmap workflow figure ---")

# Combine different heatmap types into one comprehensive figure
workflow_plots <- list()

# 1. Representative correlation heatmap
if (length(correlation_plots) > 0) {
  workflow_plots[["correlation"]] <- correlation_plots[[1]] + 
    labs(title = "A. Gene Correlation Analysis", 
         subtitle = "Co-expression patterns in representative cancer type")
}

# 2. Patient stratification heatmap
if (length(stratification_plots) > 0) {
  workflow_plots[["stratification"]] <- stratification_plots[[1]] + 
    labs(title = "B. Patient Stratification", 
         subtitle = "Clustering based on chaperone expression")
}

# 3. Pathway analysis heatmap - prefer PNG over PDF
pathway_png_file <- file.path(heatmap_paths$pathways, "pathways_heatmap.png")
pathway_pdf_file <- file.path(heatmap_paths$pathways, "pathways_heatmap.pdf")

if (file.exists(pathway_png_file)) {
  img <- read_image_for_plot(pathway_png_file)
  if (!is.null(img)) {
    workflow_plots[["pathways"]] <- create_image_plot(
      img, 
      title = "C. Pathway Enrichment",
      subtitle = "Functional pathway analysis"
    )
  }
} else if (file.exists(pathway_pdf_file)) {
  img <- read_image_for_plot(pathway_pdf_file)
  if (!is.null(img)) {
    workflow_plots[["pathways"]] <- create_image_plot(
      img, 
      title = "C. Pathway Enrichment",
      subtitle = "Functional pathway analysis"
    )
  }
}

# 4. Top enriched pathways heatmap - prefer PNG over PDF
top_pathways_png_file <- file.path(heatmap_paths$pathways, "top_enriched_pathways.png")
top_pathways_pdf_file <- file.path(heatmap_paths$pathways, "top_enriched_pathways.pdf")

if (file.exists(top_pathways_png_file)) {
  img <- read_image_for_plot(top_pathways_png_file)
  if (!is.null(img)) {
    workflow_plots[["top_pathways"]] <- create_image_plot(
      img, 
      title = "D. Top Enriched Pathways",
      subtitle = "Most significant pathway associations"
    )
  }
} else if (file.exists(top_pathways_pdf_file)) {
  img <- read_image_for_plot(top_pathways_pdf_file)
  if (!is.null(img)) {
    workflow_plots[["top_pathways"]] <- create_image_plot(
      img, 
      title = "D. Top Enriched Pathways",
      subtitle = "Most significant pathway associations"
    )
  }
}

if (length(workflow_plots) >= 3) {
  # Create a comprehensive workflow figure
  if (length(workflow_plots) == 4) {
    # 2x2 layout
    workflow_figure <- wrap_plots(workflow_plots, ncol = 2, nrow = 2)
  } else {
    # 1x3 layout
    workflow_figure <- wrap_plots(workflow_plots, ncol = 3, nrow = 1)
  }
  
  workflow_figure <- workflow_figure + 
    plot_annotation(
      title = "Comprehensive Chaperone Gene Heatmap Analysis",
      subtitle = "Multi-level heatmap analysis from correlation to pathway enrichment",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14)
      )
    )
  
  ggsave(file.path(output_dir, "Figure_3_Comprehensive_heatmap_workflow.png"), 
         workflow_figure, width = 20, height = 12, dpi = 300, bg = "white")
  
  ggsave(file.path(output_dir, "Figure_3_Comprehensive_heatmap_workflow.pdf"), 
         workflow_figure, width = 20, height = 12, bg = "white")
  
  message("Created comprehensive heatmap workflow figure")
}

# =============================================================================
# APPROACH 4: Side-by-side comparison figure
# =============================================================================

message("\n--- Creating side-by-side comparison figure ---")

# Create before/after or comparison figures
comparison_plots <- list()

# Compare different cancer types side by side
cancer_pairs <- list(
  c("TCGA-BRCA", "TCGA-LUAD"),
  c("TCGA-COAD", "TCGA-PRAD"),
  c("TCGA-HNSC", "TCGA-LIHC")
)

for (i in seq_along(cancer_pairs)) {
  pair <- cancer_pairs[[i]]
  pair_plots <- list()
  
  for (cancer in pair) {
    cancer_file <- correlation_files[grep(cancer, correlation_files)]
    if (length(cancer_file) > 0) {
      img <- read_image_for_plot(cancer_file[1])
      if (!is.null(img)) {
        plot_title <- gsub("TCGA-", "", cancer)
        pair_plots[[cancer]] <- create_image_plot(img, title = plot_title)
      }
    }
  }
  
  if (length(pair_plots) == 2) {
    comparison_plots[[paste0("pair_", i)]] <- wrap_plots(pair_plots, ncol = 2, nrow = 1)
  }
}

if (length(comparison_plots) > 0) {
  final_comparison <- wrap_plots(comparison_plots, ncol = 1, nrow = length(comparison_plots))
  
  final_comparison <- final_comparison + 
    plot_annotation(
      title = "Cancer Type Comparison: Chaperone Gene Correlation Patterns",
      subtitle = "Side-by-side comparison of correlation patterns across cancer types",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )
  
  ggsave(file.path(output_dir, "Figure_4_Cancer_type_comparison.png"), 
         final_comparison, width = 16, height = 18, dpi = 300, bg = "white")
  
  ggsave(file.path(output_dir, "Figure_4_Cancer_type_comparison.pdf"), 
         final_comparison, width = 16, height = 18, bg = "white")
  
  message("Created cancer type comparison figure")
}

# =============================================================================
# Summary and recommendations
# =============================================================================

message("\n=== Summary of Generated Heatmap Figures ===")
message("Generated figures saved in: ", output_dir)
message("Heatmap figure types created:")
message("- Figure 1: Multi-cancer correlation heatmaps (4x2 grid)")
message("- Figure 2: Patient stratification heatmaps (3x2 grid)")
message("- Figure 3: Comprehensive heatmap workflow (mixed layout)")
message("- Figure 4: Cancer type heatmap comparison (side-by-side)")

# Create a summary document
summary_file <- file.path(output_dir, "Heatmap_Publication_Figure_Guide.md")
summary_content <- "# Heatmap Publication Figure Guide

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
"

writeLines(summary_content, summary_file)

message("Heatmap publication figure guide created: ", summary_file)
message("Script completed at: ", Sys.time())

# Close log
sink()
close(con)
