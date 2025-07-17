library(readr)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(WGCNA)
library(stringr)
library(tidyr)
library(ggplot2)

# WGCNA parameters
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)

# Define file paths
rds_dir <- "TCGA-Chaperones/rds"
cluster_dir <- "TCGA-Chaperones/heatmap_patients_stratification/deseq2_clustering"
output_dir <- "TCGA-Chaperones/WGCNA_by_clusters/processed_data"
log_file <- "TCGA-Chaperones/WGCNA_by_clusters/normalization_log.txt"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
con <- file(log_file, "w")
sink(con, split = TRUE)

# Function to normalize and prepare data for WGCNA by clusters
normalize_by_clusters <- function(project_id) {
  message("\n====== Processing ", project_id, " for cluster-based WGCNA ======")
  
  # File paths
  cluster_file <- file.path(cluster_dir, project_id, paste0(project_id, "_patient_clusters.csv"))
  
  # Check if cluster file exists
  if (!file.exists(cluster_file)) {
    message("  Cluster file not found for ", project_id, ". Skipping.")
    return(NULL)
  }
  
  # Read cluster assignments
  message("  Loading cluster assignments")
  cluster_data <- read_csv(cluster_file, show_col_types = FALSE)
  
  # Ensure Patient_Cluster is character
  cluster_data$Patient_Cluster <- as.character(cluster_data$Patient_Cluster)
  
  # Get unique clusters
  unique_clusters <- sort(unique(cluster_data$Patient_Cluster))
  n_clusters <- length(unique_clusters)
  
  message("  Found ", n_clusters, " patient clusters: ", paste(unique_clusters, collapse = ", "))
  
  # Show cluster distribution
  cluster_counts <- table(cluster_data$Patient_Cluster)
  message("  Cluster distribution:")
  print(cluster_counts)
  
  # Find and load the RDS file for this project
  rds_files <- list.files(rds_dir, pattern = paste0("^", project_id, "_.*\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("  No RDS file found for project ", project_id, ". Skipping.")
    return(NULL)
  }
  
  message("  Loading RDS file: ", basename(rds_files[1]))
  se_object <- readRDS(rds_files[1])
  
  # Filter for tumor samples only
  message("  Filtering for tumor samples")
  if ("definition" %in% colnames(colData(se_object))) {
    tumor_samples <- se_object$definition == "Primary solid Tumor"
    se_object <- se_object[, tumor_samples]
    message("  Kept ", sum(tumor_samples), " tumor samples")
  }
  
  # Get count data
  assay_names <- assayNames(se_object)
  message("  Available assays: ", paste(assay_names, collapse = ", "))
  
  # Choose the appropriate assay for expression data
  if ("unstranded" %in% assay_names) {
    count_data <- assay(se_object, "unstranded")
    message("  Using unstranded counts")
  } else if ("stranded_first" %in% assay_names) {
    count_data <- assay(se_object, "stranded_first")
    message("  Using stranded_first counts")
  } else if ("stranded_second" %in% assay_names) {
    count_data <- assay(se_object, "stranded_second")
    message("  Using stranded_second counts")
  } else if ("counts" %in% assay_names) {
    count_data <- assay(se_object, "counts")
    message("  Using counts")
  } else {
    message("  No suitable count data found. Skipping.")
    return(NULL)
  }
  
  # Get patient barcodes and create short versions for matching
  tcga_barcodes <- colnames(count_data)
  short_barcodes <- substr(tcga_barcodes, 1, 15)
  
  # Create mapping for cluster data
  cluster_data$short_barcode <- substr(cluster_data$Patient_ID, 1, 15)
  
  # Process each cluster separately
  results <- list()
  
  for (cluster_id in unique_clusters) {
    message("\n  === Processing Cluster ", cluster_id, " ===")
    
    # Get patients in this cluster
    cluster_patients <- cluster_data$short_barcode[cluster_data$Patient_Cluster == cluster_id]
    
    # Find matching indices in RDS data
    cluster_indices <- which(short_barcodes %in% cluster_patients)
    
    if (length(cluster_indices) < 10) {
      message("  Not enough samples in cluster ", cluster_id, " (", length(cluster_indices), " samples). Skipping.")
      next
    }
    
    message("  Found ", length(cluster_indices), " samples in cluster ", cluster_id)
    
    # Extract count data for this cluster
    cluster_counts <- count_data[, cluster_indices]
    
    # Round counts to integers (DESeq2 requirement)
    cluster_counts <- round(cluster_counts)
    
    # Create sample metadata
    sample_metadata <- data.frame(
      sample_id = colnames(cluster_counts),
      cluster = paste0("Cluster_", cluster_id),
      row.names = colnames(cluster_counts)
    )
    
    # Create DESeq2 dataset
    message("  Creating DESeq2 dataset")
    dds <- DESeqDataSetFromMatrix(
      countData = cluster_counts,
      colData = sample_metadata,
      design = ~ 1  # Intercept only design
    )
    
    # Filter low-count genes
    # Keep genes with at least 10 counts in at least 25% of samples
    min_samples <- ceiling(0.25 * ncol(dds))
    keep <- rowSums(counts(dds) >= 10) >= min_samples
    dds <- dds[keep,]
    
    message("  Filtered ", sum(keep), " genes (from ", nrow(cluster_counts), " total)")
    
    # Run DESeq2 normalization
    message("  Running DESeq2 normalization")
    dds <- DESeq(dds)
    
    # Apply variance stabilizing transformation
    message("  Applying variance stabilizing transformation")
    vsd <- vst(dds, blind = TRUE)
    vst_data <- assay(vsd)
    
    # Transpose for WGCNA (samples as rows, genes as columns)
    wgcna_data <- t(vst_data)
    
    # Quality control: check for outliers
    message("  Performing quality control")
    sample_clustering <- hclust(dist(wgcna_data), method = "average")
    
    # Create output directory for this project
    project_output_dir <- file.path(output_dir, project_id)
    dir.create(project_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Save cluster data
    cluster_file_name <- paste0(project_id, "_cluster_", cluster_id, "_wgcna_data.rds")
    saveRDS(wgcna_data, file.path(project_output_dir, cluster_file_name))
    
    # Save sample clustering plot
    pdf(file.path(project_output_dir, paste0(project_id, "_cluster_", cluster_id, "_sample_clustering.pdf")), 
        width = 12, height = 8)
    plot(sample_clustering, 
         main = paste(project_id, "- Cluster", cluster_id, "Sample Clustering"),
         xlab = "Sample Index",
         ylab = "Height",
         cex = 0.8)
    dev.off()
    
    # Calculate expression statistics
    mean_expr <- apply(wgcna_data, 2, mean)
    var_expr <- apply(wgcna_data, 2, var)
    
    expr_summary <- data.frame(
      Gene = colnames(wgcna_data),
      Mean = mean_expr,
      Variance = var_expr,
      CV = sqrt(var_expr) / mean_expr
    )
    
    # Save expression summary
    write_csv(expr_summary, file.path(project_output_dir, paste0(project_id, "_cluster_", cluster_id, "_expression_summary.csv")))
    
    # Store results
    results[[paste0("Cluster_", cluster_id)]] <- list(
      cluster_id = cluster_id,
      n_samples = nrow(wgcna_data),
      n_genes = ncol(wgcna_data),
      mean_expression = mean(mean_expr),
      file_path = file.path(project_output_dir, cluster_file_name)
    )
    
    message("  Cluster ", cluster_id, " processed successfully")
    message("    Samples: ", nrow(wgcna_data))
    message("    Genes: ", ncol(wgcna_data))
    message("    Mean expression: ", round(mean(mean_expr), 3))
  }
  
  # Save processing summary
  if (length(results) > 0) {
    summary_df <- data.frame(
      Project = project_id,
      Cluster = sapply(results, function(x) x$cluster_id),
      N_Samples = sapply(results, function(x) x$n_samples),
      N_Genes = sapply(results, function(x) x$n_genes),
      Mean_Expression = sapply(results, function(x) x$mean_expression),
      Data_File = sapply(results, function(x) basename(x$file_path))
    )
    
    write_csv(summary_df, file.path(project_output_dir, paste0(project_id, "_processing_summary.csv")))
    message("\n  Processing summary saved for ", project_id)
    
    return(summary_df)
  } else {
    message("\n  No clusters were successfully processed for ", project_id)
    return(NULL)
  }
}

# Find all project directories in the clustering output
project_dirs <- list.dirs(cluster_dir, full.names = FALSE, recursive = FALSE)
project_dirs <- project_dirs[project_dirs != ""]

message("Found ", length(project_dirs), " projects with cluster data")
print(project_dirs)

# Process each project
all_results <- list()

for (project_id in project_dirs) {
  result <- try(normalize_by_clusters(project_id))
  if (!inherits(result, "try-error") && !is.null(result)) {
    all_results[[project_id]] <- result
  }
}

# Generate overall summary
if (length(all_results) > 0) {
  overall_summary <- bind_rows(all_results)
  write_csv(overall_summary, file.path(output_dir, "overall_processing_summary.csv"))
  
  message("\n====== Overall Summary ======")
  message("Projects processed: ", length(all_results))
  message("Total clusters processed: ", nrow(overall_summary))
  message("Average samples per cluster: ", round(mean(overall_summary$N_Samples), 1))
  message("Average genes per cluster: ", round(mean(overall_summary$N_Genes), 1))
  
  print(overall_summary)
} else {
  message("No projects were successfully processed")
}

# Close log connection
sink()
close(con)

message("\nCluster-based data normalization for WGCNA complete!")
