# Load the RDS file
combined_results <- readRDS(file.path("TCGA-WGCNA", "TCGA_gene_expression_analysis_results.rds"))
# Extract the main components
expression_data <- combined_results$expression_data
match_results <- combined_results$match_results

# Examine the structure without loading all data into view
str(combined_results, max.level = 2)

# Get the names of all datasets
dataset_names <- names(expression_data)
cat("Available datasets:", length(dataset_names), "\n")
print(dataset_names)

# Function to explore a specific dataset
explore_dataset <- function(dataset_name) {
  if (!(dataset_name %in% dataset_names)) {
    stop("Dataset not found")
  }
  
  # Get expression data for this dataset
  expr_data <- expression_data[[dataset_name]]
  
  # Get match results for this dataset
  matches <- match_results[[dataset_name]]
  
  # Print summary information
  cat("Dataset:", dataset_name, "\n")
  cat("Expression matrix dimensions:", dim(expr_data)[1], "genes x", dim(expr_data)[2], "samples\n")
  cat("Matching genes:", matches$match_count, "\n")
  
  # Return a small sample of the expression data (first 5 genes, first 5 samples)
  if (nrow(expr_data) > 0 && ncol(expr_data) > 0) {
    sample_size <- min(5, nrow(expr_data))
    sample_cols <- min(5, ncol(expr_data))
    return(expr_data[1:sample_size, 1:sample_cols])
  } else {
    return(NULL)
  }
}

# Save individual datasets as separate RDS files
save_individual_datasets <- function(output_dir = (file.path("TCGA-WGCNA", "split_datasets"))) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Save each dataset as a separate file
  for (name in dataset_names) {
    # Create a clean filename
    clean_name <- gsub("\\.rds$", "", name)
    filename <- file.path(output_dir, paste0(clean_name, "_processed.rds"))
    
    # Create a smaller object with just this dataset's information
    dataset_obj <- list(
      expression_data = expression_data[[name]],
      match_info = match_results[[name]]
    )
    
    # Save to file
    saveRDS(dataset_obj, filename)
    cat("Saved:", filename, "\n")
  }
}

# Example usage:
# 1. To see list of available datasets:
print(dataset_names)

# 2. To explore a specific dataset (replace "DATASET_NAME" with actual name):
# sample_data <- explore_dataset("DATASET_NAME")
# View(sample_data)

# 3. To save all datasets as individual files:
save_individual_datasets()