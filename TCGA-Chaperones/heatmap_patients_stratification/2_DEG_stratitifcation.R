library(readr)
library(dplyr)
library(SummarizedExperiment)
library(edgeR)
library(stringr)
library(tidyr)
library(ggplot2)

# Define file paths
rds_dir <- "TCGA-Chaperones/rds"
expr_dir <- "TCGA-Chaperones/heatmap_patients_stratification/heatmaps"
output_dir <- "TCGA-Chaperones/heatmap_patients_stratification/DEG_results"
log_file <- "TCGA-Chaperones/heatmap_patients_stratification/DEG_analysis_log.txt"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
con <- file(log_file, "w")
sink(con, split = TRUE)

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
}

# Function to perform DEG analysis for a project
analyze_project <- function(project_id) {
  message("\n====== Processing ", project_id, " ======")
  
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
  
  # Choose the appropriate assay for expression data - prioritizing counts for DEG
  if ("counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "counts")
    message("  Using counts for differential expression analysis")
  } else if ("HTSeq - Counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "HTSeq - Counts")
    message("  Using HTSeq - Counts for differential expression analysis")
  } else if ("raw_counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "raw_counts")
    message("  Using raw_counts for differential expression analysis")
  } else if ("unstranded" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "unstranded")
    message("  Using unstranded counts for differential expression analysis")
  } else {
    message("  No count data found. Cannot perform differential expression analysis.")
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
  
  # Create DGEList object
  sample_indices <- c(high_indices, low_indices)
  group <- factor(c(rep("High", length(high_indices)), 
                   rep("Low", length(low_indices))), 
                 levels = c("High", "Low"))
  
  # Create a DGEList object
  dge <- DGEList(counts = count_data[, sample_indices], 
                 group = group)
  
  # Filter low expressed genes
  message("  Filtering lowly expressed genes")
  keep <- filterByExpr(dge)
  message("  Keeping ", sum(keep), " out of ", nrow(dge), " genes")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalize the data
  message("  Normalizing data")
  dge <- calcNormFactors(dge)
  
  # Design matrix
  design <- model.matrix(~group)
  
  # Estimate dispersion
  message("  Estimating dispersion")
  dge <- estimateDisp(dge, design)
  
  # Fit the model
  message("  Fitting GLM")
  fit <- glmQLFit(dge, design)
  
  # Test for differential expression
  message("  Testing for differential expression")
  qlf <- glmQLFTest(fit, coef = 2)  # High vs. Low
  
  # Extract results
  message("  Extracting DEG results")
  res <- topTags(qlf, n = Inf)
  results_df <- as.data.frame(res)
  results_df$gene_ensembl <- rownames(results_df)
  
  # Add gene names if available in the SummarizedExperiment object
    if (is(se_object, "SummarizedExperiment") && is(try(rowRanges(se_object), silent = TRUE), "GRanges")) {
    gene_info <- as.data.frame(mcols(rowRanges(se_object)))
    
    # Find column with gene names
    potential_name_cols <- c("gene_name", "external_gene_name", "symbol", "gene_symbol")
    name_col <- potential_name_cols[potential_name_cols %in% colnames(gene_info)][1]
    
    if (!is.na(name_col)) {
      # Create mapping from ENSEMBL ID to gene name
      ensembl_to_gene_map <- data.frame(
        ensembl_id = rownames(se_object),
        gene_name = gene_info[[name_col]],
        stringsAsFactors = FALSE
      )
      
      # Extract base ENSEMBL IDs
      ensembl_to_gene_map$ensembl_base <- extract_base_ensembl(ensembl_to_gene_map$ensembl_id)
      results_df$ensembl_base <- extract_base_ensembl(results_df$gene_ensembl)
      
      # Merge gene names to results
      results_df <- merge(results_df, 
                         ensembl_to_gene_map[, c("ensembl_base", "gene_name")], 
                         by = "ensembl_base", all.x = TRUE)
    }
  }
  
  # Reorder columns to put gene names first if they exist
  if ("gene_name" %in% colnames(results_df)) {
    results_df <- results_df %>%
      select(ensembl_base, gene_name, gene_ensembl, everything()) %>%
      arrange(PValue)
  } else {
    results_df <- results_df %>%
      select(ensembl_base, gene_ensembl, everything()) %>%
      arrange(PValue)
  }
  
  # Save results
  results_file <- file.path(output_dir, paste0(project_id, "_high_vs_low_DEG.csv"))
  write_csv(results_df, results_file)
  
  # Create summary of results
  up_genes <- sum(results_df$logFC > 0 & results_df$FDR < 0.05)
  down_genes <- sum(results_df$logFC < 0 & results_df$FDR < 0.05)
  total_sig <- sum(results_df$FDR < 0.05)
  
  message("  DEG analysis complete. Results saved to ", results_file)
  message("  Summary: ", total_sig, " significant genes (FDR < 0.05)")
  message("           ", up_genes, " up-regulated in high expression group")
  message("           ", down_genes, " down-regulated in high expression group")
  
  # Create MA plot if there are significant DEGs
  if (total_sig > 0) {
    plot_file <- file.path(output_dir, paste0(project_id, "_high_vs_low_MAplot.pdf"))
    pdf(plot_file, width = 10, height = 8)
    
    # Create data frame for plotting
    plot_data <- results_df %>%
      mutate(significance = case_when(
        FDR < 0.01 ~ "FDR < 0.01",
        FDR < 0.05 ~ "FDR < 0.05",
        TRUE ~ "Not significant"
      ),
      significance = factor(significance, levels = c("FDR < 0.01", "FDR < 0.05", "Not significant")))
    
    # Create MA plot
    p <- ggplot(plot_data, aes(x = logCPM, y = logFC, color = significance)) +
      geom_point(size = 0.8, alpha = 0.7) +
      scale_color_manual(values = c("FDR < 0.01" = "red", "FDR < 0.05" = "orange", "Not significant" = "grey")) +
      labs(title = paste(project_id, "- High vs Low Chaperone Expression"),
           subtitle = paste("Up:", up_genes, "Down:", down_genes, "Total Sig:", total_sig),
           x = "Average Expression (logCPM)",
           y = "Log Fold Change (High vs Low)") +
      theme_bw() +
      theme(legend.position = "bottom")
    
    print(p)
    dev.off()
    message("  MA plot saved to ", plot_file)
    
    # Also create a volcano plot
    plot_file <- file.path(output_dir, paste0(project_id, "_high_vs_low_volcano.pdf"))
    pdf(plot_file, width = 10, height = 8)
    
    # Create volcano plot
    p <- ggplot(plot_data, aes(x = logFC, y = -log10(FDR), color = significance)) +
      geom_point(size = 0.8, alpha = 0.7) +
      scale_color_manual(values = c("FDR < 0.01" = "red", "FDR < 0.05" = "orange", "Not significant" = "grey")) +
      labs(title = paste(project_id, "- High vs Low Chaperone Expression"),
           subtitle = paste("Up:", up_genes, "Down:", down_genes, "Total Sig:", total_sig),
           x = "Log Fold Change (High vs Low)",
           y = "-log10 FDR") +
      theme_bw() +
      theme(legend.position = "bottom") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey")
    
    print(p)
    dev.off()
    message("  Volcano plot saved to ", plot_file)
  }
  
  # Return a summary of the results
  return(data.frame(
    Project = project_id,
    High_Samples = length(high_indices),
    Low_Samples = length(low_indices),
    Filtered_Genes = sum(keep),
    Total_DEGs = total_sig,
    Up_in_High = up_genes,
    Down_in_High = down_genes,
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
  result <- try(analyze_project(project_id))
  if (!inherits(result, "try-error") && !is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Generate summary table
if (length(results) > 0) {
  summary_df <- bind_rows(results)
  write_csv(summary_df, file.path(output_dir, "DEG_analysis_summary.csv"))
  
  message("\n====== Overall Summary ======")
  message("Total projects processed: ", nrow(summary_df))
  message("Total projects with DEGs: ", sum(summary_df$Total_DEGs > 0))
  message("Total significant DEGs: ", sum(summary_df$Total_DEGs))
  message("Average DEGs per project: ", round(mean(summary_df$Total_DEGs), 1))
  
  # Print summary table
  print(summary_df)
} else {
  message("No projects were successfully processed")
}

# Close log connection
sink()
close(con)

message("\nDifferential expression analysis complete!")