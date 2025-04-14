library(dplyr)
library(readr)
library(stringr)
library(SummarizedExperiment)

# Base directory for DEG results
base_dir <- "TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/DEG_results"

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  # Split by period and take the first part (base ID without version)
  return(sapply(strsplit(ensembl_id, "\\."), function(x) x[1]))
}

# Get all cancer projects (subdirectories)
projects <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
message(paste("Found", length(projects), "projects with DEG results"))

# Process all DEG result files for all projects
total_files <- 0
files_processed <- 0

# Create a mapping cache to avoid re-reading the same RDS files
ensembl_to_gene_cache <- list()

for (project in projects) {
  project_path <- file.path(base_dir, project)
  message(paste("Processing project:", project))
  
  # Get all DEG result files in this project
  deg_files <- list.files(project_path, pattern = "_DEG_results\\.csv$", full.names = TRUE)
  total_files <- total_files + length(deg_files)
  
  if (length(deg_files) == 0) {
    message(paste("  No DEG result files found in", project))
    next
  }
  
  # Load the gene name mapping for this project if we haven't already
  if (!project %in% names(ensembl_to_gene_cache)) {
    # Find the corresponding RDS file
    rds_path <- list.files("TCGA-Chaperones/rds", 
                          pattern = paste0("^", project, "_.*\\.rds$"), 
                          full.names = TRUE)
    
    if (length(rds_path) == 0) {
      message(paste("  Warning: Could not find RDS file for project", project))
      next
    }
    
    message(paste("  Loading gene metadata from", basename(rds_path[1])))
    
    # Load the RDS file to extract gene metadata
    se_object <- readRDS(rds_path[1])
    
    # Extract gene metadata from SummarizedExperiment object
    if (is(rowRanges(se_object), "GRanges")) {
      gene_info <- as.data.frame(mcols(rowRanges(se_object)))
      
      # Check which columns might contain gene names
      potential_name_cols <- c("gene_name", "external_gene_name", "symbol", "gene_symbol")
      name_col <- potential_name_cols[potential_name_cols %in% colnames(gene_info)][1]
      
      if (!is.na(name_col)) {
        # Create mapping from ENSEMBL ID to gene name
        ensembl_to_gene_map <- data.frame(
          ensembl_id = rownames(se_object),
          gene_name = gene_info[[name_col]]
        )
        
        # Extract base ENSEMBL IDs without version
        ensembl_to_gene_map$ensembl_base <- extract_base_ensembl(ensembl_to_gene_map$ensembl_id)
        
        # Store in cache
        ensembl_to_gene_cache[[project]] <- ensembl_to_gene_map
        message(paste("  Found", nrow(ensembl_to_gene_map), "gene name mappings"))
      } else {
        message("  Warning: Could not find gene name column in metadata")
        next
      }
    } else {
      message("  Warning: RDS file does not have expected gene metadata structure")
      next
    }
  }
  
  # Get the mapping for this project
  ensembl_to_gene_map <- ensembl_to_gene_cache[[project]]
  
  # Process each DEG file
  for (file_path in deg_files) {
    file_name <- basename(file_path)
    message(paste("  Processing file:", file_name))
    
    # Read the DEG results
    deg_data <- read_csv(file_path, show_col_types = FALSE)
    
    # Extract base ENSEMBL IDs (without version numbers)
    if ("gene_ensembl" %in% colnames(deg_data)) {
      deg_data$ensembl_base <- extract_base_ensembl(deg_data$gene_ensembl)
      
      # FIX: Use base R instead of dplyr select
      mapping_subset <- ensembl_to_gene_map[, c("ensembl_base", "gene_name")]
      deg_data <- merge(deg_data, mapping_subset, by="ensembl_base", all.x=TRUE)
      
      # Report how many IDs were successfully mapped
      mapped_count <- sum(!is.na(deg_data$gene_name))
      message(paste("    Successfully mapped", mapped_count, "out of", nrow(deg_data), "entries"))
      
      # Save the updated data
      write_csv(deg_data, file_path)
      files_processed <- files_processed + 1
    } else {
      message("    No 'gene_ensembl' column found in this file, skipping")
    }
  }
}

message(paste("Processing complete!", files_processed, "files processed out of", total_files, "total files"))