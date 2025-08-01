# Create Comprehensive Pan-Cancer SummarizedExperiment
# This script loads all TCGA .rds files, extracts ALL genes (not just chaperones),
# and combines them into a single SummarizedExperiment for comprehensive analysis

# Load required libraries
library(tibble)
library(readr)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(stringr)

# Define file paths
raw_data_dir <- "raw_data"
output_dir <- "TCGA-Chaperones/Patient_stratification/pan_cancer_comprehensive"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Output file for the comprehensive SummarizedExperiment
comprehensive_se_file <- file.path(output_dir, "pan_cancer_all_genes_SE.rds")

message("=== Creating Comprehensive Pan-Cancer SummarizedExperiment ===")
message("Starting at: ", Sys.time())

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
}

# Function to process a single TCGA project and extract ALL genes
process_tcga_project_all_genes <- function(project_id) {
  message("  Processing project: ", project_id)
  
  # Find and load RDS file from raw_data directory
  rds_files <- list.files(raw_data_dir, pattern = paste0("^", project_id, "_.*transcriptomic_exp\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("    No RDS file found for project ", project_id)
    return(NULL)
  }
  
  se_object <- readRDS(rds_files[1])
  
  # Filter for tumor samples only
  if ("definition" %in% colnames(colData(se_object))) {
    se_object <- se_object[, se_object$definition == "Primary solid Tumor"]
  } else if ("sample_type" %in% colnames(colData(se_object))) {
    se_object <- se_object[, se_object$sample_type == "Primary Tumor"]
  } else if (any(grepl("shortLetterCode", colnames(colData(se_object))))) {
    se_object <- se_object[, se_object$shortLetterCode == "TP"]
  }
  
  n_samples <- ncol(se_object)
  message("    Tumor samples: ", n_samples)
  
  if (n_samples < 10) {
    message("    Too few samples - skipping")
    return(NULL)
  }
  
  # Get RAW COUNT data for DESeq2
  assay_names <- assayNames(se_object)
  
  if ("unstranded" %in% assay_names) {
    count_data <- assay(se_object, "unstranded")
    used_assay <- "unstranded"
  } else if ("stranded_first" %in% assay_names) {
    count_data <- assay(se_object, "stranded_first")
    used_assay <- "stranded_first"
  } else if ("stranded_second" %in% assay_names) {
    count_data <- assay(se_object, "stranded_second")
    used_assay <- "stranded_second"
  } else {
    count_data <- assay(se_object, assay_names[1])
    used_assay <- assay_names[1]
    message("    WARNING: Using ", assay_names[1], " (may not be raw counts)")
  }
  
  # Round counts to integers
  count_data <- round(count_data)
  
  # Add project prefix to sample names
  colnames(count_data) <- paste0(project_id, ":", colnames(count_data))
  
  # Create enhanced sample metadata
  sample_metadata <- data.frame(
    Patient_ID = colnames(count_data),
    Project = project_id,
    Sample_ID = str_extract(colnames(count_data), "[^:]+$"),
    sample_type = "tumor",
    used_assay = used_assay,
    stringsAsFactors = FALSE,
    row.names = NULL  # Explicitly set to NULL to avoid warnings
  )
  rownames(sample_metadata) <- colnames(count_data)
  
  message("    Total genes: ", nrow(count_data))
  message("    Used assay: ", used_assay)
  
  return(list(
    count_data = count_data,
    sample_metadata = sample_metadata,
    project = project_id,
    n_samples = n_samples,
    n_genes = nrow(count_data),
    used_assay = used_assay
  ))
}

# Function to combine all TCGA projects with ALL genes
combine_all_tcga_projects <- function(force_recompute = FALSE) {
  
  # Check if comprehensive SE already exists
  if (!force_recompute && file.exists(comprehensive_se_file)) {
    message("\n=== Loading existing comprehensive SummarizedExperiment ===")
    comprehensive_se <- readRDS(comprehensive_se_file)
    message("Loaded comprehensive SE with ", ncol(comprehensive_se), " samples and ", nrow(comprehensive_se), " genes")
    return(comprehensive_se)
  }
  
  message("\n=== Processing all TCGA projects for comprehensive analysis ===")
  
  # Get all TCGA .rds files from raw_data directory
  tcga_files <- list.files(raw_data_dir, pattern = "^TCGA-.*_transcriptomic_exp\\.rds$", full.names = TRUE)
  message("Found ", length(tcga_files), " TCGA files")
  
  # Extract project IDs from filenames
  extract_project_id <- function(filename) {
    basename(filename) %>%
      str_replace("_transcriptomic_exp\\.rds$", "")
  }
  
  project_ids <- sapply(tcga_files, extract_project_id)
  message("Projects to process: ", paste(project_ids, collapse = ", "))
  
  # Process all projects
  all_count_data <- list()
  all_sample_metadata <- list()
  project_summaries <- list()
  
  for (i in seq_along(project_ids)) {
    project_id <- project_ids[i]
    
    tryCatch({
      project_result <- process_tcga_project_all_genes(project_id)
      if (!is.null(project_result)) {
        all_count_data[[project_id]] <- project_result$count_data
        all_sample_metadata[[project_id]] <- project_result$sample_metadata
        project_summaries[[project_id]] <- list(
          project = project_id,
          n_samples = project_result$n_samples,
          n_genes = project_result$n_genes,
          used_assay = project_result$used_assay
        )
      }
    }, error = function(e) {
      message("    Error processing ", project_id, ": ", e$message)
    })
  }
  
  if (length(all_count_data) == 0) {
    stop("No projects were successfully processed!")
  }
  
  message("\nSuccessfully processed ", length(all_count_data), " projects")
  
  # Print project summaries
  for (proj in names(project_summaries)) {
    summary <- project_summaries[[proj]]
    message("  ", proj, ": ", summary$n_samples, " samples, ", summary$n_genes, " genes (", summary$used_assay, ")")
  }
  
  # Combine all count matrices
  message("\n=== Combining count matrices ===")
  
  # Get common genes across all projects (intersection)
  all_gene_names <- Reduce(intersect, lapply(all_count_data, rownames))
  message("Common genes across all projects: ", length(all_gene_names))
  
  if (length(all_gene_names) == 0) {
    stop("No common genes found across projects!")
  }
  
  # Combine count matrices for common genes
  combined_count_matrices <- list()
  for (project_id in names(all_count_data)) {
    count_matrix <- all_count_data[[project_id]]
    combined_count_matrices[[project_id]] <- count_matrix[all_gene_names, , drop = FALSE]
  }
  
  # Merge all count matrices
  combined_count_data <- do.call(cbind, combined_count_matrices)
  
  # Combine all sample metadata
  combined_sample_metadata <- do.call(rbind, all_sample_metadata)
  
  message("Combined count matrix: ", nrow(combined_count_data), " genes × ", ncol(combined_count_data), " samples")
  message("Combined sample metadata: ", nrow(combined_sample_metadata), " samples")
  
  # Ensure row names match
  rownames(combined_sample_metadata) <- combined_sample_metadata$Patient_ID
  
  # Create comprehensive SummarizedExperiment
  message("\n=== Creating comprehensive SummarizedExperiment ===")
  
  comprehensive_se <- SummarizedExperiment(
    assays = list(counts = combined_count_data),
    colData = S4Vectors::DataFrame(combined_sample_metadata),
    metadata = list(
      creation_date = Sys.time(),
      n_projects = length(project_summaries),
      project_summaries = project_summaries,
      common_genes = all_gene_names,
      description = "Comprehensive pan-cancer TCGA SummarizedExperiment with all common genes"
    )
  )
  
  # Add gene information to rowData if possible
  gene_info <- data.frame(
    gene_id = all_gene_names,
    base_gene_id = extract_base_ensembl(all_gene_names),
    stringsAsFactors = FALSE
  )
  rownames(gene_info) <- all_gene_names
  SummarizedExperiment::rowData(comprehensive_se) <- S4Vectors::DataFrame(gene_info)
  
  message("Created comprehensive SummarizedExperiment:")
  message("  Samples: ", ncol(comprehensive_se))
  message("  Genes: ", nrow(comprehensive_se))
  message("  Projects: ", length(project_summaries))
  
  # Save comprehensive SummarizedExperiment
  message("\nSaving comprehensive SummarizedExperiment to: ", comprehensive_se_file)
  saveRDS(comprehensive_se, comprehensive_se_file)
  
  return(comprehensive_se)
}

# Main execution function
main <- function(force_recompute = FALSE) {
  
  # Create comprehensive SummarizedExperiment
  comprehensive_se <- combine_all_tcga_projects(force_recompute = force_recompute)
  
  # Print final summary
  message("\n=== Comprehensive Pan-Cancer SummarizedExperiment Created! ===")
  message("Output file: ", comprehensive_se_file)
  message("Total samples: ", ncol(comprehensive_se))
  message("Total genes: ", nrow(comprehensive_se))
  message("Projects included: ", length(S4Vectors::metadata(comprehensive_se)$project_summaries))
  
  # Show some basic statistics
  project_counts <- table(comprehensive_se$Project)
  message("\nSamples per project:")
  for (i in seq_along(project_counts)) {
    message("  ", names(project_counts)[i], ": ", project_counts[i])
  }
  
  message("\nFile ready for comprehensive DEG analysis!")
  
  return(comprehensive_se)
}

# Run the analysis
# Set force_recompute=TRUE to recreate the comprehensive SE from scratch
result <- main(force_recompute = TRUE)

message("\n=== Comprehensive SE creation completed at: ", Sys.time(), " ===")
