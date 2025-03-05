# Get list of all .rds files in the raw_data folder
rds_files <- list.files(path = "raw_data", pattern = "\\.rds$", full.names = TRUE)

# Check if any files were found
if(length(rds_files) == 0) {
  stop("No .rds files found in raw_data folder")
}

# Print the files found
cat("Found", length(rds_files), "RDS files in raw_data folder:\n")
print(rds_files)

# Load the ProteinID.csv file
protein_ids <- read.csv("TCGA-gene_expression/ProteinID.csv", stringsAsFactors = FALSE)
cat("Loaded protein IDs with", nrow(protein_ids), "entries\n")
cat("ENSEMBL IDs to match against:", paste(head(protein_ids$ENSEMBL), collapse=", "), "...\n\n")

# Create a list to store results for each dataset
results_list <- list()
# Create a list to track gene matches
match_results <- list()

# Loop through each file
for(file_path in rds_files) {
  # Extract just the filename for reference
  file_name <- basename(file_path)
  cat("\nProcessing file:", file_name, "\n")
  
  # Read the data
  data <- readRDS(file_path)
  
  # Process if the data has the expected structure
  tryCatch({
    # Convert unstranded to counts if present
    if("unstranded" %in% assayNames(data)) {
      assayNames(data)[assayNames(data) == "unstranded"] <- "counts"
    }
    
    # Check sample types
    cat("Sample types in this dataset:\n")
    print(table(data$sample_type))
    
    # Filter to keep only Primary Tumor samples
    if("Primary Tumor" %in% data$sample_type) {
      data <- data[, data$sample_type == "Primary Tumor"]
      
      # Create DESeq dataset with filtered data
      dds <- DESeqDataSet(data, design = ~ 1)
      dds <- DESeq(dds)
      vsd <- varianceStabilizingTransformation(dds)
      
      # Continue with analysis
      wpn_vsd <- getVarianceStabilizedData(dds)
      rv_wpn <- rowVars(wpn_vsd)
      
      q75_wpn <- quantile(rowVars(wpn_vsd), .75) 
      q95_wpn <- quantile(rowVars(wpn_vsd), .95)
      expr_normalized <- wpn_vsd[rv_wpn > q95_wpn, ]
      
      # Check for matches between gene names in expr_normalized and ENSEMBL column
      gene_names <- rownames(expr_normalized)
      
      # Clean ENSEMBL IDs to remove version numbers if present (e.g., ENSG00000123456.7 -> ENSG00000123456)
      clean_genes <- sub("\\.\\d+$", "", gene_names)
      clean_ensembl <- sub("\\.\\d+$", "", protein_ids$ENSEMBL)
      
      # Find matches
      matches <- intersect(clean_genes, clean_ensembl)
      match_count <- length(matches)
      
      # Print match results
      cat("\nNumber of matches with protein list:", match_count, "\n")
      if (match_count > 0) {
        cat("Matching genes:", paste(head(matches, 5), collapse=", "))
        if (match_count > 5) cat(", ...")
        cat("\n")
      }
      
      # Store match results
      match_results[[file_name]] <- list(
        has_matches = match_count > 0,
        match_count = match_count,
        matches = matches
      )
      
      # Store results
      results_list[[file_name]] <- expr_normalized
    } else {
      cat("No 'Primary Tumor' samples found in this dataset. Skipping...\n")
      match_results[[file_name]] <- list(has_matches = FALSE, match_count = 0, matches = character(0))
    }
  }, error = function(e) {
    cat("Error processing file", file_name, ":", conditionMessage(e), "\n")
    match_results[[file_name]] <- list(has_matches = FALSE, match_count = 0, matches = character(0))
  })
}

# Print summary of projects with and without matches
cat("\n\n=== GENE MATCH SUMMARY ===\n")
cat("Projects with matches:\n")
with_matches <- names(match_results)[sapply(match_results, function(x) x$has_matches)]
if (length(with_matches) > 0) {
  for (proj in with_matches) {
    cat("  - ", proj, ": ", match_results[[proj]]$match_count, " matches\n", sep="")
  }
} else {
  cat("  None\n")
}

cat("\nProjects without matches:\n")
without_matches <- names(match_results)[!sapply(match_results, function(x) x$has_matches)]
if (length(without_matches) > 0) {
  cat("  ", paste(without_matches, collapse=", "), "\n")
} else {
  cat("  None\n")
}

# Combine results into a single object
combined_results <- list(
  expression_data = results_list,
  match_results = match_results
)

# Save the combined results to an RDS file
output_file <- file.path("TCGA-WGCNA", "TCGA_gene_expression_analysis_results.rds")
saveRDS(combined_results, file = output_file)
cat("\nResults saved to:", output_file, "\n")