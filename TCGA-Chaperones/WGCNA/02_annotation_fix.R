elibrary(readr)
library(dplyr)
library(SummarizedExperiment)

# Define file paths
rds_dir <- "TCGA-Chaperones/rds"
wgcna_dir <- "TCGA-Chaperones/WGCNA/processed_data"
output_dir <- "TCGA-Chaperones/WGCNA/gene_annotation"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Find processed projects from WGCNA
wgcna_files <- list.files(wgcna_dir, pattern = "_all_expr_wgcna.rds$", full.names = FALSE)
project_ids <- unique(gsub("_all_expr_wgcna.rds$", "", wgcna_files))

message("Found ", length(project_ids), " projects to process")

# Function to extract gene names from RDS files
extract_gene_names <- function(project_id) {
  message("\n====== Processing ", project_id, " ======")
  
  # Find the RDS file
  rds_files <- list.files(rds_dir, pattern = paste0("^", project_id, "_.*\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("  No RDS file found for project ", project_id)
    return(NULL)
  }
  
  message("  Loading RDS file: ", basename(rds_files[1]))
  se_object <- readRDS(rds_files[1])
  
  # Get gene IDs from the SummarizedExperiment object
  gene_ids <- rownames(se_object)
  gene_symbols <- rep(NA, length(gene_ids))
  gene_names <- rep(NA, length(gene_ids))
  
  message("  Found ", length(gene_ids), " genes in the dataset")
  
  # Extract base Ensembl IDs (without version)
  base_ids <- gsub("\\.[0-9]+$", "", gene_ids)
  
  # Try to find gene information in different places
  if (is(rowData(se_object), "DataFrame") && ncol(rowData(se_object)) > 0) {
    message("  Extracting gene information from rowData")
    gene_info <- as.data.frame(rowData(se_object))
    
    # Print the column names to help with debugging
    message("  Available columns in rowData: ", paste(colnames(gene_info), collapse=", "))
    
    # Check for gene symbols in common column names
    symbol_cols <- c("gene_name", "symbol", "external_gene_name", "gene_symbol", "hgnc_symbol",
                    "Symbol", "gene", "GeneName")
    for (col in symbol_cols) {
      if (col %in% colnames(gene_info)) {
        gene_symbols <- gene_info[[col]]
        message("  Found gene symbols in column: ", col)
        break
      }
    }
    
    # Check for gene names/descriptions in common column names
    name_cols <- c("gene_description", "description", "full_name", "name", 
                  "gene_full_name", "gene_desc", "Description")
    for (col in name_cols) {
      if (col %in% colnames(gene_info)) {
        gene_names <- gene_info[[col]]
        message("  Found gene descriptions in column: ", col)
        break
      }
    }
  } else {
    message("  No rowData information found")
  }
  
  # Create a data frame with the information we've found
  gene_mapping <- data.frame(
    ensembl_id = gene_ids,
    ensembl_base_id = base_ids,
    gene_symbol = gene_symbols,
    gene_description = gene_names,
    stringsAsFactors = FALSE
  )
  
  # Check if we found any symbols
  symbol_count <- sum(!is.na(gene_symbols))
  message("  Found symbols for ", symbol_count, " out of ", length(gene_ids), " genes (", 
          round(symbol_count/length(gene_ids)*100, 1), "%)")
  
  # Load the WGCNA expression file to get genes actually used in analysis
  wgcna_file <- file.path(wgcna_dir, paste0(project_id, "_all_expr_wgcna.rds"))
  if (file.exists(wgcna_file)) {
    message("  Loading WGCNA preprocessed data to match genes")
    wgcna_data <- readRDS(wgcna_file)
    wgcna_genes <- colnames(wgcna_data)
    
    # Filter gene mapping to only include genes used in WGCNA
    gene_mapping_filtered <- gene_mapping %>%
      filter(ensembl_id %in% wgcna_genes | ensembl_base_id %in% gsub("\\.[0-9]+$", "", wgcna_genes))
    
    message("  Filtered to ", nrow(gene_mapping_filtered), " genes used in WGCNA")
    
    # Save the filtered mapping
    output_file_filtered <- file.path(output_dir, paste0(project_id, "_wgcna_gene_mapping.csv"))
    write_csv(gene_mapping_filtered, output_file_filtered)
    message("  Saved filtered gene mapping to: ", basename(output_file_filtered))
  }
  
  # Save the full mapping
  output_file <- file.path(output_dir, paste0(project_id, "_full_gene_mapping.csv"))
  write_csv(gene_mapping, output_file)
  message("  Saved complete gene mapping to: ", basename(output_file))
  
  return(data.frame(
    Project = project_id,
    Total_Genes = length(gene_ids),
    Genes_With_Symbols = symbol_count,
    Symbol_Coverage_Percent = round(symbol_count/length(gene_ids)*100, 1),
    stringsAsFactors = FALSE
  ))
}

# Process each project
results <- list()
for (project_id in project_ids) {
  result <- try(extract_gene_names(project_id))
  if (!inherits(result, "try-error") && !is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Generate summary table
if (length(results) > 0) {
  summary_df <- bind_rows(results)
  write_csv(summary_df, file.path(output_dir, "gene_mapping_summary.csv"))
  
  message("\n====== Overall Summary ======")
  message("Total projects processed: ", nrow(summary_df))
  message("Average gene symbol coverage: ", round(mean(summary_df$Symbol_Coverage_Percent), 1), "%")
  
  # Print summary table
  print(summary_df)
} else {
  message("No projects were successfully processed")
}

message("\nGene name extraction complete!")