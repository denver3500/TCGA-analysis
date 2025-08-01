# Comprehensive DEG Analysis: Cluster 4 vs Cluster 7 - ALL GENES
# This script loads the comprehensive SummarizedExperiment and performs
# differential expression analysis on ALL genes using patient IDs from clustering

# Load required libraries
library(tibble)
library(readr)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(stringr)
library(limma)

# Define file paths
comprehensive_se_file <- "TCGA-Chaperones/Patient_stratification/pan_cancer_comprehensive/pan_cancer_all_genes_SE.rds"
pictures_dir <- "TCGA-Chaperones/Patient_stratification/pictures"
output_dir <- "TCGA-Chaperones/Patient_stratification/comprehensive_DEG_cluster_4_vs_7"
tcga_reference_file <- "raw_data/TCGA-ACC_transcriptomic_exp.rds"  # Reference for gene names

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("=== Comprehensive DEG Analysis: Cluster 4 vs Cluster 7 (ALL GENES) ===")
message("Starting at: ", Sys.time())

# Function to load gene name mapping from TCGA reference file
load_gene_name_mapping <- function() {
  message("\n=== Loading gene name mapping ===")
  
  if (!file.exists(tcga_reference_file)) {
    stop("TCGA reference file not found: ", tcga_reference_file)
  }
  
  # Load reference TCGA file to get gene metadata
  tcga_ref <- readRDS(tcga_reference_file)
  gene_metadata <- rowData(tcga_ref)
  
  # Create mapping from Ensembl ID to gene name
  gene_mapping <- data.frame(
    ensembl_id = gene_metadata$gene_id,
    gene_name = gene_metadata$gene_name,
    stringsAsFactors = FALSE
  )
  
  # Remove rows with missing gene names
  gene_mapping <- gene_mapping[!is.na(gene_mapping$gene_name) & gene_mapping$gene_name != "", ]
  
  message("Loaded gene name mapping for ", nrow(gene_mapping), " genes")
  
  return(gene_mapping)
}

# Function to load comprehensive SummarizedExperiment
load_comprehensive_se <- function() {
  message("\n=== Loading comprehensive SummarizedExperiment ===")
  
  if (!file.exists(comprehensive_se_file)) {
    stop("Comprehensive SummarizedExperiment not found. Please run create_comprehensive_SE.R first.")
  }
  
  comprehensive_se <- readRDS(comprehensive_se_file)
  
  message("Loaded comprehensive SE:")
  message("  Samples: ", ncol(comprehensive_se))
  message("  Genes: ", nrow(comprehensive_se))
  message("  Projects: ", length(S4Vectors::metadata(comprehensive_se)$project_summaries))
  
  return(comprehensive_se)
}

# Function to load cluster patient IDs
load_cluster_patient_ids <- function() {
  message("\n=== Loading cluster patient IDs ===")
  
  cluster_4_file <- file.path(pictures_dir, "cluster_4_patient_ids.txt")
  cluster_7_file <- file.path(pictures_dir, "cluster_7_patient_ids.txt")
  
  if (!file.exists(cluster_4_file) || !file.exists(cluster_7_file)) {
    stop("Cluster patient ID files not found. Please run pan_cancer_stratification.R first.")
  }
  
  cluster_4_patients <- readLines(cluster_4_file)
  cluster_7_patients <- readLines(cluster_7_file)
  
  message("Cluster 4 patients: ", length(cluster_4_patients))
  message("Cluster 7 patients: ", length(cluster_7_patients))
  
  return(list(
    cluster_4 = cluster_4_patients,
    cluster_7 = cluster_7_patients
  ))
}

# Function to extract cluster samples from comprehensive SE
extract_cluster_samples <- function(comprehensive_se, patient_clusters) {
  message("\n=== Extracting cluster samples from comprehensive SE ===")
  
  # Get available patient IDs in the SE
  available_patients <- colnames(comprehensive_se)
  message("Available patients in comprehensive SE: ", length(available_patients))
  
  # Find matching patients for each cluster
  cluster_4_found <- patient_clusters$cluster_4[patient_clusters$cluster_4 %in% available_patients]
  cluster_7_found <- patient_clusters$cluster_7[patient_clusters$cluster_7 %in% available_patients]
  
  message("Cluster 4 patients found in SE: ", length(cluster_4_found))
  message("Cluster 7 patients found in SE: ", length(cluster_7_found))
  
  if (length(cluster_4_found) == 0 || length(cluster_7_found) == 0) {
    stop("No patients found for one or both clusters in the comprehensive SE!")
  }
  
  # Extract subset of SE for these patients
  cluster_patients <- c(cluster_4_found, cluster_7_found)
  cluster_se <- comprehensive_se[, cluster_patients]
  
  # Add cluster information to colData
  cluster_se$Cluster <- factor(
    c(rep("Cluster_4", length(cluster_4_found)), 
      rep("Cluster_7", length(cluster_7_found))),
    levels = c("Cluster_4", "Cluster_7")
  )
  
  message("Created cluster SE with ", ncol(cluster_se), " samples and ", nrow(cluster_se), " genes")
  
  return(list(
    cluster_se = cluster_se,
    cluster_4_patients = cluster_4_found,
    cluster_7_patients = cluster_7_found
  ))
}

# Function to perform comprehensive DEG analysis using DESeq2
perform_comprehensive_deg <- function(cluster_data, gene_mapping) {
  message("\n=== Performing comprehensive DEG analysis with DESeq2 ===")
  
  cluster_se <- cluster_data$cluster_se
  
  # Filter low count genes (more stringent for genome-wide analysis)
  message("Filtering low count genes...")
  min_samples <- ceiling(0.1 * ncol(cluster_se))  # At least 10% of samples
  min_count <- 10
  
  keep <- rowSums(assay(cluster_se, "counts") >= min_count) >= min_samples
  cluster_se_filtered <- cluster_se[keep, ]
  
  message("Genes after filtering: ", nrow(cluster_se_filtered), " (removed ", 
          nrow(cluster_se) - nrow(cluster_se_filtered), " low-count genes)")
  
  # Create DESeq2 dataset
  message("Creating DESeq2 dataset...")
  dds <- DESeqDataSet(cluster_se_filtered, design = ~ Cluster)
  
  # Run DESeq2 analysis
  message("Running DESeq2 analysis...")
  dds <- DESeq(dds, quiet = FALSE)
  
  # Get results
  message("Extracting results...")
  res <- results(dds, contrast = c("Cluster", "Cluster_7", "Cluster_4"))
  
  # Convert to data frame and add gene information
  deg_results <- as.data.frame(res)
  deg_results$Ensembl_ID <- rownames(deg_results)
  
  # Add gene names by merging with gene mapping
  deg_results <- merge(gene_mapping, deg_results, by.x = "ensembl_id", by.y = "Ensembl_ID", all.y = TRUE)
  
  # Replace missing gene names with Ensembl IDs
  deg_results$gene_name[is.na(deg_results$gene_name)] <- deg_results$ensembl_id[is.na(deg_results$gene_name)]
  
  # Reorder columns: Gene_Name, Ensembl_ID, then other columns
  deg_results <- deg_results[, c("gene_name", "ensembl_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  colnames(deg_results)[1:2] <- c("Gene_Name", "Ensembl_ID")
  
  # Add significance categories
  deg_results$Significance <- "Not Significant"
  deg_results$Significance[!is.na(deg_results$padj) & deg_results$padj < 0.05 & deg_results$log2FoldChange > 0.5] <- "Up in Cluster 7"
  deg_results$Significance[!is.na(deg_results$padj) & deg_results$padj < 0.05 & deg_results$log2FoldChange < -0.5] <- "Up in Cluster 4"
  
  # Remove genes with NA adjusted p-values
  deg_results <- deg_results[!is.na(deg_results$padj), ]
  
  # Sort by adjusted p-value
  deg_results <- deg_results[order(deg_results$padj), ]
  
  message("DEG analysis complete:")
  message("  Total genes analyzed: ", nrow(deg_results))
  message("  Significantly up in Cluster 7: ", sum(deg_results$Significance == "Up in Cluster 7"))
  message("  Significantly up in Cluster 4: ", sum(deg_results$Significance == "Up in Cluster 4"))
  message("  Not significant: ", sum(deg_results$Significance == "Not Significant"))
  
  # Save results
  deg_file <- file.path(output_dir, "comprehensive_DEG_cluster_4_vs_7_results.csv")
  write.csv(deg_results, deg_file, row.names = FALSE)
  message("Saved comprehensive DEG results to: ", deg_file)
  
  # Save significant genes separately
  significant_genes <- deg_results[deg_results$Significance != "Not Significant", ]
  if (nrow(significant_genes) > 0) {
    sig_file <- file.path(output_dir, "significant_DEG_cluster_4_vs_7.csv")
    write.csv(significant_genes, sig_file, row.names = FALSE)
    message("Saved ", nrow(significant_genes), " significant genes to: ", sig_file)
  }
  
  return(list(
    deg_results = deg_results,
    dds = dds,
    significant_genes = significant_genes
  ))
}

# Function to generate enhanced volcano plot
generate_enhanced_volcano <- function(deg_results) {
  message("\n=== Generating enhanced volcano plot ===")
  
  volcano_file <- file.path(output_dir, "comprehensive_volcano_plot.pdf")
  pdf(volcano_file, width = 12, height = 10)
  
  tryCatch({
    # Prepare data for plotting
    plot_data <- deg_results
    plot_data$neg_log10_padj <- -log10(plot_data$padj)
    
    # Set colors
    colors <- ifelse(plot_data$Significance == "Up in Cluster 7", "red",
                    ifelse(plot_data$Significance == "Up in Cluster 4", "blue", "grey"))
    
    # Create the plot
    plot(plot_data$log2FoldChange, plot_data$neg_log10_padj,
         col = colors, pch = 16, cex = 0.6,
         xlab = "Log2 Fold Change (Cluster 7 vs Cluster 4)",
         ylab = "-Log10 Adjusted P-value",
         main = paste0("Comprehensive Volcano Plot: Cluster 4 vs Cluster 7\n",
                      nrow(plot_data), " genes analyzed"))
    
    # Add significance lines
    abline(h = -log10(0.05), col = "black", lty = 2, lwd = 2)
    abline(v = c(-0.5, 0.5), col = "black", lty = 2, lwd = 2)
    
    # Add legend
    legend("topright", 
           legend = c(paste("Up in Cluster 7 (", sum(plot_data$Significance == "Up in Cluster 7"), ")"),
                     paste("Up in Cluster 4 (", sum(plot_data$Significance == "Up in Cluster 4"), ")"),
                     "Not Significant"),
           col = c("red", "blue", "grey"),
           pch = 16, cex = 0.8)
    
    # Label top 10 significant genes
    top_genes <- head(plot_data[plot_data$padj < 0.001 & abs(plot_data$log2FoldChange) > 1, ], 10)
    if (nrow(top_genes) > 0) {
      text(top_genes$log2FoldChange, top_genes$neg_log10_padj, 
           labels = top_genes$Gene_Name, pos = 3, cex = 0.6, col = "black")
    }
    
  }, error = function(e) {
    message("Error generating volcano plot: ", e$message)
    plot.new()
    text(0.5, 0.5, paste("Error:", e$message), cex = 1.2)
  })
  
  dev.off()
  message("Generated enhanced volcano plot: ", volcano_file)
}

# Function to generate top significant genes heatmap
generate_top_significant_heatmap <- function(cluster_data, deg_analysis, n_genes = 100) {
  message("\n=== Generating top significant genes heatmap ===")
  
  significant_genes <- deg_analysis$significant_genes
  
  if (nrow(significant_genes) == 0) {
    message("No significantly differentially expressed genes found!")
    return()
  }
  
  # Select top genes by adjusted p-value
  top_genes <- head(significant_genes, n_genes)
  message("Selected top ", nrow(top_genes), " significant genes for heatmap")
  
  # Get normalized counts from DESeq2
  dds <- deg_analysis$dds
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Extract expression data for top genes
  top_genes_expression <- normalized_counts[top_genes$Ensembl_ID, , drop = FALSE]
  
  # Set row names to gene names for better visualization
  rownames(top_genes_expression) <- top_genes$Gene_Name
  
  # Log2 transform and center
  top_genes_expression_log <- log2(top_genes_expression + 1)
  top_genes_expression_scaled <- t(scale(t(top_genes_expression_log)))
  
  # Create annotation
  patient_annotation <- data.frame(
    Cluster = cluster_data$cluster_se$Cluster
  )
  rownames(patient_annotation) <- colnames(top_genes_expression_scaled)
  
  # Create colors
  cluster_colors <- setNames(c("#FF6B6B", "#4ECDC4"), c("Cluster_4", "Cluster_7"))
  annotation_colors <- list(Cluster = cluster_colors)
  
  # Create color scale for scaled data
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Generate heatmap
  heatmap_file <- file.path(output_dir, "top_significant_genes_heatmap.pdf")
  pdf(heatmap_file, width = 16, height = 12)
  
  # Calculate gaps between clusters
  n_cluster_4 <- length(cluster_data$cluster_4_patients)
  cluster_gaps <- n_cluster_4
  
  tryCatch({
    pheatmap(
      top_genes_expression_scaled,
      color = color_palette,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_col = patient_annotation,
      annotation_colors = annotation_colors,
      gaps_col = cluster_gaps,
      border_color = NA,
      cellwidth = 0.05,
      cellheight = 4,
      show_rownames = TRUE,
      show_colnames = FALSE,
      fontsize_row = 4,
      fontsize = 6,
      main = paste0("Top ", nrow(top_genes), " Significant Genes\n",
                   "Cluster 4 vs Cluster 7 (Scaled Log2 Normalized Counts)"),
      filename = NA,
      na_col = "grey"
    )
  }, error = function(e) {
    message("Error generating heatmap: ", e$message)
    plot.new()
    text(0.5, 0.5, paste("Error:", e$message), cex = 1.2)
  })
  
  dev.off()
  message("Generated top significant genes heatmap: ", heatmap_file)
}

# Function to generate pathway enrichment summary
generate_deg_summary <- function(deg_analysis) {
  message("\n=== Generating DEG summary ===")
  
  deg_results <- deg_analysis$deg_results
  significant_genes <- deg_analysis$significant_genes
  
  # Create summary statistics
  summary_stats <- data.frame(
    Metric = c("Total genes analyzed",
               "Genes with valid adjusted p-values",
               "Significantly up in Cluster 7",
               "Significantly up in Cluster 4",
               "Total significant genes",
               "Most significant up in Cluster 7",
               "Most significant up in Cluster 4"),
    Value = c(nrow(deg_results),
              sum(!is.na(deg_results$padj)),
              sum(deg_results$Significance == "Up in Cluster 7"),
              sum(deg_results$Significance == "Up in Cluster 4"),
              sum(deg_results$Significance != "Not Significant"),
              ifelse(sum(deg_results$Significance == "Up in Cluster 7") > 0,
                     deg_results$Gene_Name[deg_results$Significance == "Up in Cluster 7"][1], "None"),
              ifelse(sum(deg_results$Significance == "Up in Cluster 4") > 0,
                     deg_results$Gene_Name[deg_results$Significance == "Up in Cluster 4"][1], "None")),
    stringsAsFactors = FALSE
  )
  
  # Save summary
  summary_file <- file.path(output_dir, "DEG_analysis_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  message("Saved DEG analysis summary to: ", summary_file)
  
  return(summary_stats)
}

# Main execution function
main <- function() {
  
  # Load comprehensive SummarizedExperiment
  comprehensive_se <- load_comprehensive_se()
  
  # Load gene name mapping
  gene_mapping <- load_gene_name_mapping()
  
  # Load cluster patient IDs
  patient_clusters <- load_cluster_patient_ids()
  
  # Extract cluster samples
  cluster_data <- extract_cluster_samples(comprehensive_se, patient_clusters)
  
  # Perform comprehensive DEG analysis
  deg_analysis <- perform_comprehensive_deg(cluster_data, gene_mapping)
  
  # Generate enhanced volcano plot
  generate_enhanced_volcano(deg_analysis$deg_results)
  
  # Generate top significant genes heatmap
  generate_top_significant_heatmap(cluster_data, deg_analysis)
  
  # Generate summary
  summary_stats <- generate_deg_summary(deg_analysis)
  
  # Print final summary
  message("\n=== Comprehensive DEG Analysis Complete! ===")
  message("Output directory: ", output_dir)
  message("Generated files:")
  message("- comprehensive_DEG_cluster_4_vs_7_results.csv (all genes)")
  message("- significant_DEG_cluster_4_vs_7.csv (significant genes only)")
  message("- comprehensive_volcano_plot.pdf")
  message("- top_significant_genes_heatmap.pdf")
  message("- DEG_analysis_summary.csv")
  
  print(summary_stats)
  
  return(list(
    cluster_data = cluster_data,
    deg_analysis = deg_analysis,
    summary_stats = summary_stats
  ))
}

# Run the comprehensive analysis
result <- main()

message("\n=== Comprehensive DEG analysis completed at: ", Sys.time(), " ===")
