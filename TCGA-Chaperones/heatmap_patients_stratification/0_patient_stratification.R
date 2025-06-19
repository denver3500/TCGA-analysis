# Load required libraries
library(readr)
library(dplyr)
library(SummarizedExperiment)
library(pheatmap)
library(RColorBrewer)

# Define file paths
project_counts_file <- "TCGA-Chaperones/higher_in_tumor_gene_counts.csv"
gene_list_file <- "TCGA-Chaperones/gene_list.csv"
rds_dir <- "TCGA-Chaperones/rds"
output_dir <- "TCGA-Chaperones/heatmap_patients_stratification/heatmaps"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load project data and filter for projects with >= 10 higher_in_tumor_count
projects_data <- read_csv(project_counts_file, show_col_types = FALSE) %>%
  filter(higher_in_tumor_count >= 10) 

message("Selected projects with higher_in_tumor_count >= 10:")
print(projects_data)

# Load gene list
gene_list <- read_csv(gene_list_file, show_col_types = FALSE)
message("Loaded gene list with ", nrow(gene_list), " genes")

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
}

# Function to generate heatmap for a given scaling method
generate_heatmap <- function(expr_data, project_id, scaling_method, scaling_suffix) {
  heatmap_file <- file.path(output_dir, paste0(project_id, "_chaperone_heatmap_", scaling_suffix, ".pdf"))
  
  n_patients <- ncol(expr_data)
  show_col_names <- n_patients <= 60
  
  pdf(heatmap_file, width = 16, height = 10)
  
  pheatmap(expr_data,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_method = "ward.D2",
           show_colnames = show_col_names,
           fontsize_row = 8,
           fontsize_col = 6,
           main = paste0(project_id, " - Chaperone Expression (", scaling_method, ")"),
           filename = NA)
  
  dev.off()
  
  message("  Generated heatmap (", scaling_method, "): ", heatmap_file)
  return(heatmap_file)
}

# Main function to process each project
process_project <- function(project_id) {
  message("\n=== Processing project: ", project_id, " ===")
  
  # Find and load RDS file
  rds_files <- list.files(rds_dir, pattern = paste0("^", project_id, "_.*\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("  No RDS file found for project ", project_id)
    return(NULL)
  }
  
  message("  Loading RDS file: ", basename(rds_files[1]))
  se_object <- readRDS(rds_files[1])
  
  # Filter for tumor samples only
  if ("definition" %in% colnames(colData(se_object))) {
    se_object <- se_object[, se_object$definition == "Primary solid Tumor"]
  } else if ("sample_type" %in% colnames(colData(se_object))) {
    se_object <- se_object[, se_object$sample_type == "Primary Tumor"]
  } else if (any(grepl("shortLetterCode", colnames(colData(se_object))))) {
    se_object <- se_object[, se_object$shortLetterCode == "TP"]
  }
  
  message("  Using ", ncol(se_object), " tumor samples")
  
  # Get expression data - Prioritize tpm_unstrand
  assay_names <- assayNames(se_object)
  message("  Available assays: ", paste(assay_names, collapse = ", "))
  
  used_assay <- NULL
  if ("tpm_unstrand" %in% assay_names) {
    expr_data <- assay(se_object, "tpm_unstrand")
    used_assay <- "tpm_unstrand"
    message("  Using tpm_unstrand data")
  } else if ("tpm" %in% assay_names) {
    expr_data <- assay(se_object, "tpm")
    used_assay <- "tpm"
    message("  Using tpm data")
  } else if ("fpkm_unstrand" %in% assay_names) {
    expr_data <- assay(se_object, "fpkm_unstrand")
    used_assay <- "fpkm_unstrand"
    message("  Using fpkm_unstrand data")
  } else if ("fpkm" %in% assay_names) {
    expr_data <- assay(se_object, "fpkm")
    used_assay <- "fpkm"
    message("  Using fpkm data")
  } else if ("normalized_count" %in% assay_names) {
    expr_data <- assay(se_object, "normalized_count")
    used_assay <- "normalized_count"
    message("  Using normalized_count data")
  } else {
    expr_data <- assay(se_object, assay_names[1])
    used_assay <- assay_names[1]
    message("  Using ", assay_names[1], " data")
  }
  
  # Get all gene IDs and create mapping
  all_genes <- rownames(expr_data)
  
  if (any(grepl("\\.", all_genes))) {
    # Remove version numbers from ENSEMBL IDs
    all_base_genes <- extract_base_ensembl(all_genes)
    gene_mapping <- data.frame(
      full_id = all_genes,
      base_id = all_base_genes,
      stringsAsFactors = FALSE
    )
  } else {
    gene_mapping <- data.frame(
      full_id = all_genes,
      base_id = all_genes,
      stringsAsFactors = FALSE
    )
  }
  
  # Find chaperone genes in the dataset
  chaperone_genes <- c()
  chaperone_names <- c()
  
  for (i in 1:nrow(gene_list)) {
    ensembl_id <- gene_list$ENSEMBL[i]
    gene_name <- gene_list$Name[i]
    
    # Find matching gene in dataset
    match_idx <- which(gene_mapping$base_id == ensembl_id)
    
    if (length(match_idx) > 0) {
      full_id <- gene_mapping$full_id[match_idx[1]]
      chaperone_genes <- c(chaperone_genes, full_id)
      chaperone_names <- c(chaperone_names, gene_name)
    }
  }
  
  if (length(chaperone_genes) == 0) {
    message("  No chaperone genes found in dataset")
    return(NULL)
  }
  
  message("  Found ", length(chaperone_genes), " of ", nrow(gene_list), " chaperone genes")
  
  # Extract chaperone expression data
  chaperone_expr <- expr_data[chaperone_genes, , drop = FALSE]
  
  # Data transformation with debugging
  message("  Original data range: ", round(min(chaperone_expr, na.rm = TRUE), 2), " to ", round(max(chaperone_expr, na.rm = TRUE), 2))
  message("  Original data median: ", round(median(chaperone_expr, na.rm = TRUE), 2))
  
  # Log-transform if using TPM/FPKM data
  if (grepl("tpm|fpkm", tolower(used_assay))) {
    message("  Applying log2 transformation for TPM/FPKM data")
    chaperone_expr <- log2(chaperone_expr + 1)
    message("  Log-transformed range: ", round(min(chaperone_expr, na.rm = TRUE), 2), " to ", round(max(chaperone_expr, na.rm = TRUE), 2))
  } else if (max(chaperone_expr, na.rm = TRUE) > 50) {
    message("  Applying log2 transformation for high-value data")
    chaperone_expr <- log2(chaperone_expr + 1)
    message("  Log-transformed range: ", round(min(chaperone_expr, na.rm = TRUE), 2), " to ", round(max(chaperone_expr, na.rm = TRUE), 2))
  } else {
    message("  No log transformation applied")
  }
  
  # Set gene names as row names
  rownames(chaperone_expr) <- chaperone_names
  
  # Remove genes with all NA values before scaling
  chaperone_expr <- chaperone_expr[!apply(is.na(chaperone_expr), 1, all), ]
  
  message("  Final matrix before scaling: ", nrow(chaperone_expr), " genes x ", ncol(chaperone_expr), " patients")
  
  # OPTION 1: Scale genes within each patient (compare genes within patient)
  message("\n  === SCALING OPTION 1: Patient-wise scaling ===")
  chaperone_expr_scaled_patients <- scale(chaperone_expr)  # Scale columns (patients)
  
  # Handle cases where all values are the same (SD = 0)
  chaperone_expr_scaled_patients[is.na(chaperone_expr_scaled_patients)] <- 0
  
  message("    Patient-scaled range: ", round(min(chaperone_expr_scaled_patients, na.rm = TRUE), 2), " to ", round(max(chaperone_expr_scaled_patients, na.rm = TRUE), 2))
  message("    Patient-scaled mean per gene: ", round(mean(rowMeans(chaperone_expr_scaled_patients, na.rm = TRUE)), 3))
  message("    Patient-scaled mean per patient: ", round(mean(colMeans(chaperone_expr_scaled_patients, na.rm = TRUE)), 3))
  
  # OPTION 2: Scale patients within each gene (compare patients for each gene)
  message("\n  === SCALING OPTION 2: Gene-wise scaling ===")
  chaperone_expr_scaled_genes <- t(scale(t(chaperone_expr)))  # Scale rows (genes)
  
  # Handle cases where all values are the same (SD = 0)
  chaperone_expr_scaled_genes[is.na(chaperone_expr_scaled_genes)] <- 0
  
  message("    Gene-scaled range: ", round(min(chaperone_expr_scaled_genes, na.rm = TRUE), 2), " to ", round(max(chaperone_expr_scaled_genes, na.rm = TRUE), 2))
  message("    Gene-scaled mean per gene: ", round(mean(rowMeans(chaperone_expr_scaled_genes, na.rm = TRUE)), 3))
  message("    Gene-scaled mean per patient: ", round(mean(colMeans(chaperone_expr_scaled_genes, na.rm = TRUE)), 3))
  
  # OPTION 3: TRUE Global scaling (center and scale using overall dataset statistics)
  message("\n  === SCALING OPTION 3: TRUE Global scaling ===")
  overall_mean <- mean(chaperone_expr, na.rm = TRUE)
  overall_sd <- sd(as.vector(chaperone_expr), na.rm = TRUE)
  chaperone_expr_scaled_global <- (chaperone_expr - overall_mean) / overall_sd
  
  # Handle cases where SD = 0
  chaperone_expr_scaled_global[is.na(chaperone_expr_scaled_global)] <- 0
  
  message("    Global-scaled range: ", round(min(chaperone_expr_scaled_global, na.rm = TRUE), 2), " to ", round(max(chaperone_expr_scaled_global, na.rm = TRUE), 2))
  message("    Global-scaled overall mean: ", round(mean(chaperone_expr_scaled_global, na.rm = TRUE), 3))
  message("    Global-scaled overall SD: ", round(sd(as.vector(chaperone_expr_scaled_global), na.rm = TRUE), 3))
  
  # Generate all 3 heatmaps
  message("\n  === Generating 3 heatmaps ===")
  
  heatmap1 <- generate_heatmap(chaperone_expr_scaled_patients, project_id, 
                              "Patient-wise Scaling", "patient_scaled")
  
  heatmap2 <- generate_heatmap(chaperone_expr_scaled_genes, project_id, 
                              "Gene-wise Scaling", "gene_scaled")
  
  heatmap3 <- generate_heatmap(chaperone_expr_scaled_global, project_id, 
                              "Global Scaling", "global_scaled")
  
  # Calculate mean chaperone expression per patient for each scaling method
  patient_mean_expr_patients <- colMeans(chaperone_expr_scaled_patients, na.rm = TRUE)
  patient_mean_expr_genes <- colMeans(chaperone_expr_scaled_genes, na.rm = TRUE)
  patient_mean_expr_global <- colMeans(chaperone_expr_scaled_global, na.rm = TRUE)
  
  # Save patient expression data for all scaling methods
  patient_file <- file.path(output_dir, paste0(project_id, "_patient_expression_all_scalings.csv"))
  patient_df <- data.frame(
    Patient_ID = names(patient_mean_expr_patients),
    Mean_Expression_Patient_Scaled = patient_mean_expr_patients,
    Mean_Expression_Gene_Scaled = patient_mean_expr_genes,
    Mean_Expression_Global_Scaled = patient_mean_expr_global,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Mean_Expression_Gene_Scaled))
  
  write_csv(patient_df, patient_file)
  message("  Saved patient expression data for all scalings: ", patient_file)
  
  # Create summary comparing all scaling methods
  scaling_comparison <- data.frame(
    Project = project_id,
    Scaling_Method = c("Patient_Scaled", "Gene_Scaled", "Global_Scaled"),
    n_patients = rep(ncol(chaperone_expr), 3),
    n_genes = rep(nrow(chaperone_expr), 3),
    assay_used = rep(used_assay, 3),
    min_value = c(min(chaperone_expr_scaled_patients, na.rm = TRUE),
                  min(chaperone_expr_scaled_genes, na.rm = TRUE),
                  min(chaperone_expr_scaled_global, na.rm = TRUE)),
    max_value = c(max(chaperone_expr_scaled_patients, na.rm = TRUE),
                  max(chaperone_expr_scaled_genes, na.rm = TRUE),
                  max(chaperone_expr_scaled_global, na.rm = TRUE)),
    mean_value = c(mean(chaperone_expr_scaled_patients, na.rm = TRUE),
                   mean(chaperone_expr_scaled_genes, na.rm = TRUE),
                   mean(chaperone_expr_scaled_global, na.rm = TRUE)),
    sd_value = c(sd(chaperone_expr_scaled_patients, na.rm = TRUE),
                 sd(chaperone_expr_scaled_genes, na.rm = TRUE),
                 sd(chaperone_expr_scaled_global, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  # Save scaling comparison
  comparison_file <- file.path(output_dir, paste0(project_id, "_scaling_comparison.csv"))
  write_csv(scaling_comparison, comparison_file)
  message("  Saved scaling comparison: ", comparison_file)
  
  message("  Scaling comparison summary:")
  print(scaling_comparison)
  
  return(data.frame(
    Project = project_id,
    Genes_Found = length(chaperone_genes),
    Samples = ncol(chaperone_expr),
    Assay_Used = used_assay,
    stringsAsFactors = FALSE
  ))
}

# Process all projects
message("\n=== Starting processing of all projects ===")

results <- list()
all_scaling_comparisons <- list()

for (i in 1:nrow(projects_data)) {
  project_id <- projects_data$project_id[i]
  result <- process_project(project_id)
  if (!is.null(result)) {
    results[[length(results) + 1]] <- result
    
    # Also collect scaling comparison for combined file
    comparison_file <- file.path(output_dir, paste0(project_id, "_scaling_comparison.csv"))
    if (file.exists(comparison_file)) {
      scaling_comparison <- read_csv(comparison_file, show_col_types = FALSE)
      all_scaling_comparisons[[length(all_scaling_comparisons) + 1]] <- scaling_comparison
    }
  }
}

# Save overall summary
if (length(results) > 0) {
  summary_df <- bind_rows(results)
  summary_file <- file.path(output_dir, "processing_summary.csv")
  write_csv(summary_df, summary_file)
  
  message("\n=== Processing Summary ===")
  print(summary_df)
  message("\nSummary saved to: ", summary_file)
  
  # Combine all scaling comparisons into one file
  if (length(all_scaling_comparisons) > 0) {
    combined_scaling_comparison <- bind_rows(all_scaling_comparisons)
    combined_comparison_file <- file.path(output_dir, "all_projects_scaling_comparison.csv")
    write_csv(combined_scaling_comparison, combined_comparison_file)
    message("Combined scaling comparison saved to: ", combined_comparison_file)
  }
} else {
  message("\nNo projects were successfully processed")
}

message("\n=== Processing complete! ===")
message("\nFor each project, you now have 3 heatmaps:")
message("  1. *_patient_scaled.pdf - Scale patients across genes (genes comparable within patient)")
message("  2. *_gene_scaled.pdf - Scale genes across patients (patients comparable within gene)")  
message("  3. *_global_scaled.pdf - TRUE global scaling (everything comparable to overall mean)")