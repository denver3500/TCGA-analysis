# Comprehensive DEG Analysis: Cluster 2 vs Cluster 7 - ALL GENES
# This script is a copy of DEG.R but compares cluster 2 vs cluster 7 only.

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
output_dir <- "TCGA-Chaperones/Patient_stratification/comprehensive_DEG_cluster_2_vs_7"
tcga_reference_file <- "raw_data/TCGA-ACC_transcriptomic_exp.rds"  # Reference for gene names

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("=== Comprehensive DEG Analysis: Cluster 2 vs Cluster 7 (ALL GENES) ===")
message("Starting at: ", Sys.time())

# Function to load gene name mapping from TCGA reference file
load_gene_name_mapping <- function() {
  if (!file.exists(tcga_reference_file)) {
    stop("TCGA reference file not found: ", tcga_reference_file)
  }
  tcga_ref <- readRDS(tcga_reference_file)
  gene_metadata <- rowData(tcga_ref)
  gene_mapping <- data.frame(
    ensembl_id = gene_metadata$gene_id,
    gene_name = gene_metadata$gene_name,
    stringsAsFactors = FALSE
  )
  gene_mapping <- gene_mapping[!is.na(gene_mapping$gene_name) & gene_mapping$gene_name != "", ]
  return(gene_mapping)
}

# Function to load comprehensive SummarizedExperiment
load_comprehensive_se <- function() {
  if (!file.exists(comprehensive_se_file)) {
    stop("Comprehensive SummarizedExperiment not found. Please run create_comprehensive_SE.R first.")
  }
  comprehensive_se <- readRDS(comprehensive_se_file)
  return(comprehensive_se)
}

# Function to load cluster patient IDs
load_cluster_patient_ids <- function() {
  cluster_2_file <- file.path(pictures_dir, "cluster_2_patient_ids.txt")
  cluster_7_file <- file.path(pictures_dir, "cluster_7_patient_ids.txt")
  if (!file.exists(cluster_2_file) || !file.exists(cluster_7_file)) {
    stop("Cluster patient ID files not found. Please run pan_cancer_stratification.R first.")
  }
  cluster_2_patients <- readLines(cluster_2_file)
  cluster_7_patients <- readLines(cluster_7_file)
  return(list(
    cluster_2 = cluster_2_patients,
    cluster_7 = cluster_7_patients
  ))
}

# Function to extract cluster samples from comprehensive SE
extract_cluster_samples <- function(comprehensive_se, patient_clusters) {
  available_patients <- colnames(comprehensive_se)
  cluster_2_found <- patient_clusters$cluster_2[patient_clusters$cluster_2 %in% available_patients]
  cluster_7_found <- patient_clusters$cluster_7[patient_clusters$cluster_7 %in% available_patients]
  if (length(cluster_2_found) == 0 || length(cluster_7_found) == 0) {
    stop("No patients found for one or both clusters in the comprehensive SE!")
  }
  cluster_patients <- c(cluster_2_found, cluster_7_found)
  cluster_se <- comprehensive_se[, cluster_patients]
  cluster_se$Cluster <- factor(
    c(rep("Cluster_2", length(cluster_2_found)), 
      rep("Cluster_7", length(cluster_7_found))),
    levels = c("Cluster_2", "Cluster_7")
  )
  return(list(
    cluster_se = cluster_se,
    cluster_2_patients = cluster_2_found,
    cluster_7_patients = cluster_7_found
  ))
}

# Function to perform comprehensive DEG analysis using DESeq2
perform_comprehensive_deg <- function(cluster_data, gene_mapping) {
  cluster_se <- cluster_data$cluster_se
  min_samples <- ceiling(0.1 * ncol(cluster_se))
  min_count <- 10
  keep <- rowSums(assay(cluster_se, "counts") >= min_count) >= min_samples
  cluster_se_filtered <- cluster_se[keep, ]
  dds <- DESeqDataSet(cluster_se_filtered, design = ~ Cluster)
  dds <- DESeq(dds, quiet = FALSE)
  res <- results(dds, contrast = c("Cluster", "Cluster_7", "Cluster_2"))
  deg_results <- as.data.frame(res)
  deg_results$Ensembl_ID <- rownames(deg_results)
  deg_results <- merge(gene_mapping, deg_results, by.x = "ensembl_id", by.y = "Ensembl_ID", all.y = TRUE)
  deg_results$gene_name[is.na(deg_results$gene_name)] <- deg_results$ensembl_id[is.na(deg_results$gene_name)]
  deg_results <- deg_results[, c("gene_name", "ensembl_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  colnames(deg_results)[1:2] <- c("Gene_Name", "Ensembl_ID")
  deg_results$Significance <- "Not Significant"
  deg_results$Significance[!is.na(deg_results$padj) & deg_results$padj < 0.05 & deg_results$log2FoldChange > 0.5] <- "Up in Cluster 7"
  deg_results$Significance[!is.na(deg_results$padj) & deg_results$padj < 0.05 & deg_results$log2FoldChange < -0.5] <- "Up in Cluster 2"
  deg_results <- deg_results[!is.na(deg_results$padj), ]
  deg_results <- deg_results[order(deg_results$padj), ]
  deg_file <- file.path(output_dir, "comprehensive_DEG_cluster_2_vs_7_results.csv")
  write.csv(deg_results, deg_file, row.names = FALSE)
  significant_genes <- deg_results[deg_results$Significance != "Not Significant", ]
  if (nrow(significant_genes) > 0) {
    sig_file <- file.path(output_dir, "significant_DEG_cluster_2_vs_7.csv")
    write.csv(significant_genes, sig_file, row.names = FALSE)
  }
  return(list(
    deg_results = deg_results,
    dds = dds,
    significant_genes = significant_genes
  ))
}

# Main execution function
main <- function() {
  comprehensive_se <- load_comprehensive_se()
  gene_mapping <- load_gene_name_mapping()
  patient_clusters <- load_cluster_patient_ids()
  cluster_data <- extract_cluster_samples(comprehensive_se, patient_clusters)
  deg_analysis <- perform_comprehensive_deg(cluster_data, gene_mapping)
}

result <- main()
message("\n=== Comprehensive DEG analysis for Cluster 2 vs 7 completed at: ", Sys.time(), " ===")
