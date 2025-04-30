library(dplyr)
library(readr)
library(stringr)
library(SummarizedExperiment)

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  # Split by period and take the first part (base ID without version)
  return(sapply(strsplit(ensembl_id, "\\."), function(x) x[1]))
}

# Input and output file paths
input_file <- "TCGA-BRCA_BRCA1_DEG_results_debug.csv"
output_file <- "TCGA-BRCA_BRCA1_DEG_results_with_genes.csv"

# Project ID (extracted from filename)
project_id <- "TCGA-BRCA"

# Read the DEG results
message("Reading DEG file: ", input_file)
deg_data <- read_csv(input_file, show_col_types = FALSE)

# Find the corresponding RDS file
rds_dir <- "TCGA-Chaperones/rds"
rds_pattern <- paste0("^", project_id, "_.*\\.rds$")
rds_files <- list.files(rds_dir, pattern = rds_pattern, full.names = TRUE)

if (length(rds_files) == 0) {
  stop("Could not find RDS file for project ", project_id, " in directory ", rds_dir)
}

# Use the first matching RDS file
rds_path <- rds_files[1]
message("Loading gene metadata from ", basename(rds_path))

# Load the RDS file to extract gene metadata
se_object <- readRDS(rds_path)

# Extract gene metadata from SummarizedExperiment object
if (is(rowRanges(se_object), "GRanges")) {
  gene_info <- as.data.frame(mcols(rowRanges(se_object)))
  
  # Check which columns might contain gene names
  potential_name_cols <- c("gene_name", "external_gene_name", "symbol", "gene_symbol")
  name_col <- potential_name_cols[potential_name_cols %in% colnames(gene_info)][1]
  
  if (is.na(name_col)) {
    stop("Could not find gene name column in metadata")
  }
  
  # Create mapping from ENSEMBL ID to gene name
  ensembl_to_gene_map <- data.frame(
    ensembl_id = rownames(se_object),
    gene_name = gene_info[[name_col]]
  )
  
  # Extract base ENSEMBL IDs without version
  ensembl_to_gene_map$ensembl_base <- extract_base_ensembl(ensembl_to_gene_map$ensembl_id)
  message("Found ", nrow(ensembl_to_gene_map), " gene name mappings")
} else {
  stop("RDS file does not have expected gene metadata structure")
}

# Extract base ENSEMBL IDs from the DEG data
if ("gene_ensembl" %in% colnames(deg_data)) {
  deg_data$ensembl_base <- extract_base_ensembl(deg_data$gene_ensembl)
  
  # Use base R merge instead of dplyr
  mapping_subset <- ensembl_to_gene_map[, c("ensembl_base", "gene_name")]
  deg_data <- merge(deg_data, mapping_subset, by="ensembl_base", all.x=TRUE)
  
  # Report how many IDs were successfully mapped
  mapped_count <- sum(!is.na(deg_data$gene_name))
  message("Successfully mapped ", mapped_count, " out of ", nrow(deg_data), " entries")
  
  # Save the updated data
  write_csv(deg_data, output_file)
  message("Updated file saved to: ", output_file)
} else {
  stop("No 'gene_ensembl' column found in the input file")
}