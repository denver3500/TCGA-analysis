library(readr)
library(dplyr)
library(SummarizedExperiment)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)

# Define file paths
project_counts_file <- "TCGA-Chaperones/higher_in_tumor_gene_counts.csv"
gene_list_file <- "TCGA-Chaperones/gene_list.csv"
rds_dir <- "TCGA-Chaperones/rds"
output_dir <- "TCGA-Chaperones/heatmap_patients_stratification/heatmaps"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load project data and filter for projects with > 10 higher_in_tumor_count
projects_data <- read_csv(project_counts_file, show_col_types = FALSE) %>%
  filter(higher_in_tumor_count >= 10) 

message("Selected projects with higher_in_tumor_count > 10:")
print(projects_data)

# Load gene list
gene_list <- read_csv(gene_list_file, show_col_types = FALSE)
message("Loaded gene list with ", nrow(gene_list), " genes")

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
}

# Function to process a project and generate a heatmap
process_project <- function(project_id) {
  message("\nProcessing project: ", project_id)
  
  # Find the RDS file for this project
  rds_files <- list.files(rds_dir, pattern = paste0("^", project_id, "_.*\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("  No RDS file found for project ", project_id)
    return(NULL)
  }
  
  # Load the RDS file
  message("  Loading RDS file: ", basename(rds_files[1]))
  se_object <- readRDS(rds_files[1])
  
  # Filter for tumor samples
  if ("definition" %in% colnames(colData(se_object))) {
    se_object <- se_object[, se_object$definition == "Primary solid Tumor"]
  } else if ("sample_type" %in% colnames(colData(se_object))) {
    se_object <- se_object[, se_object$sample_type == "Primary Tumor"]
  } else if (any(grepl("shortLetterCode", colnames(colData(se_object))))) {
    se_object <- se_object[, se_object$shortLetterCode == "TP"]
  }
  
  message("  Using ", ncol(se_object), " tumor samples")
  
  # Get available assay names
  assay_names <- assayNames(se_object)
  message("  Available assays: ", paste(assay_names, collapse=", "))
  
  # Choose the appropriate assay for expression data
  if ("tpm" %in% assay_names) {
    expr_data <- assay(se_object, "tpm")
  } else if ("fpkm" %in% assay_names) {
    expr_data <- assay(se_object, "fpkm")
  } else if ("normalized_count" %in% assay_names) {
    expr_data <- assay(se_object, "normalized_count")
  } else if ("counts" %in% assay_names) {
    expr_data <- assay(se_object, "counts")
  } else {
    expr_data <- assay(se_object, assay_names[1])
    message("  Using ", assay_names[1], " as expression data")
  }
  
  # Get row names (gene IDs) from the expression data
  all_genes <- rownames(expr_data)
  
  # Extract base ENSEMBL IDs if they have version numbers
  if (any(grepl("\\.", all_genes))) {
    all_base_genes <- extract_base_ensembl(all_genes)
    gene_id_mapping <- data.frame(
      full_id = all_genes,
      base_id = all_base_genes,
      stringsAsFactors = FALSE
    )
  } else {
    gene_id_mapping <- data.frame(
      full_id = all_genes,
      base_id = all_genes,
      stringsAsFactors = FALSE
    )
  }
  
  # Find our genes of interest in the dataset
  our_genes <- c()
  our_gene_names <- c()
  
  for (i in 1:nrow(gene_list)) {
    ensembl_id <- gene_list$ENSEMBL[i]
    gene_name <- gene_list$Name[i]
    
    # Find matching genes in the data
    matching_rows <- which(gene_id_mapping$base_id == ensembl_id)
    
    if (length(matching_rows) > 0) {
      full_id <- gene_id_mapping$full_id[matching_rows[1]]
      our_genes <- c(our_genes, full_id)
      our_gene_names <- c(our_gene_names, gene_name)
    }
  }
  
  if (length(our_genes) == 0) {
    message("  No matching genes found in the dataset")
    return(NULL)
  }
  
  message("  Found ", length(our_genes), " of ", nrow(gene_list), " genes in the dataset")
  
  # Extract expression data for our genes of interest
  gene_expr <- expr_data[our_genes, , drop = FALSE]
  
  # Log-transform if not already on log scale
  if (median(gene_expr, na.rm = TRUE) > 50) {
    message("  Data appears to be raw counts, applying log2 transformation")
    gene_expr <- log2(gene_expr + 1)  # Add 1 to avoid log(0)
  }
  
  # Calculate mean expression per patient
  patient_mean_expr <- colMeans(gene_expr, na.rm = TRUE)
  
  # Order patients by mean expression (descending)
  ordered_patients <- names(sort(patient_mean_expr, decreasing = TRUE))
  
  # Reorder expression matrix
  gene_expr_ordered <- gene_expr[, ordered_patients, drop = FALSE]
  
  # Rename the rows to gene names
  rownames(gene_expr_ordered) <- our_gene_names
  
  # Scale the data by row (z-score)
  gene_expr_scaled <- t(scale(t(gene_expr_ordered)))
  
  # Create annotation for patients (top 25%, middle 50%, bottom 25%)
  n_patients <- length(ordered_patients)
  quartiles <- cut(1:n_patients, 
                  breaks = c(0, n_patients * 0.25, n_patients * 0.75, n_patients), 
                  labels = c("High", "Medium", "Low"))
  
  patient_anno <- data.frame(
    Expression_Group = quartiles
  )
  rownames(patient_anno) <- ordered_patients
  
  # Define colors for the annotation
  anno_colors <- list(
    Expression_Group = c(High = "#B2182B", Medium = "#F4F4F4", Low = "#2166AC")
  )
  
  # Generate the heatmap
  heatmap_file <- file.path(output_dir, paste0(project_id, "_chaperone_expression_heatmap.pdf"))
  pdf(heatmap_file, width = 12, height = 8)
  
  # Control the number of patient labels based on sample size
  show_col_names <- n_patients <= 30
  
  pheatmap(gene_expr_scaled,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_colnames = show_col_names,
           annotation_col = patient_anno,
           annotation_colors = anno_colors,
           fontsize_row = 10,
           main = paste0(project_id, " Chaperone Expression"),
           filename = NA)  # NA means don't save, we'll save with pdf()
  dev.off()
  
  message("  Generated heatmap: ", heatmap_file)
  
  # Also save the expression data for potential future use
  expr_file <- file.path(output_dir, paste0(project_id, "_chaperone_expression_data.csv"))
  expr_df <- as.data.frame(gene_expr_ordered)
  expr_df$Gene_Name <- rownames(expr_df)
  expr_df$ENSEMBL_ID <- our_genes
  expr_df <- expr_df %>% select(Gene_Name, ENSEMBL_ID, everything())
  write_csv(expr_df, expr_file)
  
  return(data.frame(
    Project = project_id,
    Genes_Found = length(our_genes),
    Samples = n_patients
  ))
}

# Process each project
results <- list()
for (i in 1:nrow(projects_data)) {
  project_id <- projects_data$project_id[i]
  result <- process_project(project_id)
  if (!is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Combine results and save summary
if (length(results) > 0) {
  summary_df <- bind_rows(results)
  write_csv(summary_df, file.path(output_dir, "project_processing_summary.csv"))
  message("\nSummary of processed projects:")
  print(summary_df)
} else {
  message("\nNo projects were successfully processed")
}

message("\nProcessing complete!")