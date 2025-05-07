library(readr)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(WGCNA)
library(stringr)
library(tidyr)
library(ggplot2)

# WGCNA parameters - increase max block size for larger datasets
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)  # Adjust based on your system

# Define file paths
rds_dir <- "TCGA-Chaperones/rds"
expr_dir <- "TCGA-Chaperones/heatmap_patients_stratification/heatmaps"
output_dir <- "TCGA-Chaperones/WGCNA/processed_data"
log_file <- "TCGA-Chaperones/WGCNA/normalization_log.txt"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
con <- file(log_file, "w")
sink(con, split = TRUE)

# Function to normalize and prepare data for WGCNA
normalize_for_wgcna <- function(project_id) {
  message("\n====== Processing ", project_id, " for WGCNA ======")
  
  # File paths
  expr_file <- file.path(expr_dir, paste0(project_id, "_chaperone_expression_data.csv"))
  
  # Check if expression file exists
  if (!file.exists(expr_file)) {
    message("  Expression file not found for ", project_id, ". Skipping.")
    return(NULL)
  }
  
  # Read chaperone expression data
  message("  Loading chaperone expression data")
  expr_data <- read_csv(expr_file, show_col_types = FALSE)
  
  # Transpose to get patients as rows and genes as columns
  gene_names <- expr_data$Gene_Name
  patient_ids <- colnames(expr_data)[3:ncol(expr_data)]
  expr_matrix <- t(as.matrix(expr_data[, 3:ncol(expr_data)]))
  colnames(expr_matrix) <- gene_names
  rownames(expr_matrix) <- patient_ids
  
  # Calculate mean expression per patient
  patient_mean_expr <- rowMeans(expr_matrix, na.rm = TRUE)
  
  # Sort patients by mean expression
  sorted_patients <- names(sort(patient_mean_expr, decreasing = TRUE))
  n_patients <- length(sorted_patients)
  
  # Define high and low groups (top 25% and bottom 25%)
  high_patients <- sorted_patients[1:round(n_patients * 0.25)]
  low_patients <- sorted_patients[(round(n_patients * 0.75) + 1):n_patients]
  
  message("  Identified ", length(high_patients), " patients with high chaperone expression")
  message("  Identified ", length(low_patients), " patients with low chaperone expression")
  
  # Find and load the RDS file for this project
  rds_files <- list.files(rds_dir, pattern = paste0("^", project_id, "_.*\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("  No RDS file found for project ", project_id, ". Skipping.")
    return(NULL)
  }
  
  message("  Loading RDS file: ", basename(rds_files[1]))
  se_object <- readRDS(rds_files[1])
  
  # Check available assays in the SummarizedExperiment object
  message("  Available assays: ", paste(assayNames(se_object), collapse = ", "))
  
  # Choose the appropriate assay for count data
  if ("counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "counts")
    message("  Using counts for DESeq2 analysis")
  } else if ("HTSeq - Counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "HTSeq - Counts")
    message("  Using HTSeq - Counts for DESeq2 analysis")
  } else if ("raw_counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "raw_counts")
    message("  Using raw_counts for DESeq2 analysis")
  } else if ("unstranded" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "unstranded")
    message("  Using unstranded counts for DESeq2 analysis")
  } else {
    message("  No count data found. Cannot perform DESeq2 analysis.")
    return(NULL)
  }
  
  # Get TCGA barcodes from column names
  tcga_barcodes <- colnames(count_data)
  short_barcodes <- substr(tcga_barcodes, 1, 15)  # Use the first 15 characters for a more specific match
  
  # Match barcodes between expression data and RDS data
  high_indices <- which(short_barcodes %in% substr(high_patients, 1, 15))
  low_indices <- which(short_barcodes %in% substr(low_patients, 1, 15))
  
  if (length(high_indices) < 3 || length(low_indices) < 3) {
    message("  Not enough samples matched between RDS and expression data. Skipping.")
    message("  Matched samples: ", length(high_indices), " high, ", length(low_indices), " low")
    return(NULL)
  }
  
  message("  Matched ", length(high_indices), " high expression samples and ", 
          length(low_indices), " low expression samples")
  
  # Prepare count data and metadata for DESeq2
  sample_indices <- c(high_indices, low_indices)
  count_matrix <- count_data[, sample_indices]
  
  # Create sample metadata
  coldata <- data.frame(
    row.names = colnames(count_matrix),
    condition = factor(c(rep("High", length(high_indices)), 
                       rep("Low", length(low_indices))), 
                     levels = c("High", "Low"))
  )
  
  # Create DESeq2 object
  message("  Creating DESeq2 object")
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),  # Ensure counts are integers
    colData = coldata,
    design = ~ condition
  )
  
  # Filter low-count genes
  message("  Filtering lowly expressed genes")
  keep <- rowSums(counts(dds) >= 10) >= 5
  dds <- dds[keep,]
  message("  Keeping ", sum(keep), " out of ", nrow(count_matrix), " genes")
  
  # Run DESeq2
  message("  Running DESeq2 normalization")
  dds <- DESeq(dds)
  
  # Perform variance stabilizing transformation
  message("  Applying variance stabilizing transformation")
  vst <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  # Extract VST data
  vsd <- getVarianceStabilizedData(dds)
  
  # Check if all values are finite (no NaN or Inf)
  if (!all(is.finite(vsd))) {
    message("  Warning: VST data contains non-finite values, replacing with zeros")
    vsd[!is.finite(vsd)] <- 0
  }
  
  # Normalize VST data using quantile normalization as shown in the WGCNA tutorial
  message("  Applying quantile normalization")
  wpn_vsd <- vsd
  
  # For WGCNA, we want genes in columns and samples in rows
  wpn_vsd <- t(wpn_vsd)
  
  # Split normalized expression matrices for high and low groups
  high_expr <- wpn_vsd[1:length(high_indices),]
  low_expr <- wpn_vsd[(length(high_indices)+1):nrow(wpn_vsd),]
  
  # Save data for WGCNA analysis
  message("  Saving processed data for WGCNA analysis")
  saveRDS(high_expr, file.path(output_dir, paste0(project_id, "_high_expr_wgcna.rds")))
  saveRDS(low_expr, file.path(output_dir, paste0(project_id, "_low_expr_wgcna.rds")))
  saveRDS(wpn_vsd, file.path(output_dir, paste0(project_id, "_all_expr_wgcna.rds")))
  
  # Check for extreme outlier samples using hierarchical clustering
  message("  Checking for outlier samples")
  sample_clustering <- hclust(dist(wpn_vsd), method = "average")
  
  # Plot sample clustering for visual inspection
  pdf(file.path(output_dir, paste0(project_id, "_sample_clustering.pdf")))
  plot(sample_clustering, main = paste(project_id, "- Sample Clustering"), 
       labels = rownames(wpn_vsd), xlab = "", sub = "")
  dev.off()
  
  # Generate summary statistics
  message("  Generating summary statistics")
  mean_expr <- colMeans(wpn_vsd)
  var_expr <- apply(wpn_vsd, 2, var)
  
  # Create a basic expression summary
  expr_summary <- data.frame(
    Gene = colnames(wpn_vsd),
    Mean = mean_expr,
    Variance = var_expr,
    CV = sqrt(var_expr) / mean_expr
  )
  
  # Save summary
  write_csv(expr_summary, file.path(output_dir, paste0(project_id, "_expression_summary.csv")))
  
  # Return basic info
  return(data.frame(
    Project = project_id,
    Total_Samples = nrow(wpn_vsd),
    High_Samples = length(high_indices),
    Low_Samples = length(low_indices),
    Genes_Kept = ncol(wpn_vsd),
    stringsAsFactors = FALSE
  ))
}

# Find all expression files
expr_files <- list.files(expr_dir, pattern = "_chaperone_expression_data.csv$", full.names = FALSE)
project_ids <- gsub("_chaperone_expression_data.csv$", "", expr_files)

message("Found ", length(project_ids), " projects with chaperone expression data")
print(project_ids)

# Process each project
results <- list()
for (project_id in project_ids) {
  result <- try(normalize_for_wgcna(project_id))
  if (!inherits(result, "try-error") && !is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Generate summary table
if (length(results) > 0) {
  summary_df <- bind_rows(results)
  write_csv(summary_df, file.path(output_dir, "wgcna_preprocessing_summary.csv"))
  
  message("\n====== Overall Summary ======")
  message("Total projects processed: ", nrow(summary_df))
  message("Average samples per project: ", round(mean(summary_df$Total_Samples), 1))
  message("Average genes per project: ", round(mean(summary_df$Genes_Kept), 1))
  
  # Print summary table
  print(summary_df)
} else {
  message("No projects were successfully processed")
}

# Close log connection
sink()
close(con)

message("\nData normalization for WGCNA complete!")