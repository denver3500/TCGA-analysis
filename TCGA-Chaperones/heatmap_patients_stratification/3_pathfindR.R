# Load required packages
library(pathfindR)
library(dplyr)
library(readr)
library(ggkegg)
library(ggplot2)

# Define input and output directories
deg_dir <- "TCGA-Chaperones/heatmap_patients_stratification/DEG_results"
output_dir <- "TCGA-Chaperones/heatmap_patients_stratification/pathfindR_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Find all DEG files
deg_files <- list.files(deg_dir, pattern = "_high_vs_low_DEG\\.csv$", full.names = TRUE)
project_ids <- gsub(".*/(TCGA-[A-Z]+)_high_vs_low_DEG\\.csv$", "\\1", deg_files)

cat("Found", length(deg_files), "DEG files to process\n")

# Create a summary dataframe
summary_results <- list()

# Process each DEG file
for (i in 1:length(deg_files)) {
  project_id <- project_ids[i]
  deg_file <- deg_files[i]
  
  cat("\n===== Processing", project_id, "=====\n")
  start_time <- Sys.time()
  
  # Create project folder
  project_dir <- file.path(output_dir, project_id)
  dir.create(project_dir, showWarnings = FALSE)
  
  # Read and process DEG data
  deg_data <- read_csv(deg_file, show_col_types = FALSE)
  cat("Reading file:", basename(deg_file), "\n")
  
  # Filter for significant DEGs
  filtered_data <- deg_data %>% filter(FDR < 0.05)
  sig_count <- nrow(filtered_data)
  cat("Found", sig_count, "significant DEGs\n")
  
  # Skip if not enough DEGs
  if (sig_count < 10) {
    cat("Not enough significant DEGs (minimum 10). Skipping.\n")
    summary_results[[i]] <- list(
      Project = project_id,
      Significant_DEGs = sig_count,
      Enriched_Terms = 0,
      Top_Term = "Insufficient DEGs",
      Top_Term_PValue = NA,
      Runtime_Minutes = 0
    )
    next
  }
  
  # Create input for pathfindR
  input_df <- data.frame(
    Gene_symbol = filtered_data$gene_name,
    logFC = filtered_data$logFC,
    FDR_p = filtered_data$FDR
  )
  
  # Run pathfindR
  cat("Running pathfindR analysis...\n")
  res <- tryCatch({
    pathfindR_results <- run_pathfindR(
      input_df,
      gene_sets = "GO-All",
      pin_name = "KEGG",
      min_gset_size = 10,
      max_gset_size = 300,
      p_val_threshold = 0.05,
      output_dir = project_dir
    )
    
    # Save results
    result_file <- file.path(project_dir, paste0(project_id, "_pathfindR_results.csv"))
    write_csv(pathfindR_results, result_file)
    
    # Get term count and top term
    term_count <- nrow(pathfindR_results)
    if (term_count > 0) {
      top_term <- pathfindR_results$Term_Description[1]
      p_value <- pathfindR_results$lowest_p[1]
      cat("Found", term_count, "enriched terms\n")
      cat("Top term:", top_term, "(p =", p_value, ")\n")
      
      # Generate visualizations - only the ones requested
      cat("Clustering enriched terms...\n")
      cluster_enriched_terms(pathfindR_results, output_dir = project_dir)
      
      list(
        Project = project_id,
        Significant_DEGs = sig_count,
        Enriched_Terms = term_count,
        Top_Term = top_term,
        Top_Term_PValue = p_value,
        Runtime_Minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
      )
    } else {
      cat("No enriched terms found\n")
      list(
        Project = project_id,
        Significant_DEGs = sig_count, 
        Enriched_Terms = 0,
        Top_Term = "No enriched terms",
        Top_Term_PValue = NA,
        Runtime_Minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
      )
    }
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    list(
      Project = project_id,
      Significant_DEGs = sig_count,
      Enriched_Terms = 0,
      Top_Term = paste("Error:", substr(e$message, 1, 50)),
      Top_Term_PValue = NA,
      Runtime_Minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
    )
  })
  
  summary_results[[i]] <- res
}

# Convert results to dataframe
summary_df <- do.call(rbind.data.frame, lapply(summary_results, as.data.frame))

# Save summary
summary_file <- file.path(output_dir, "pathfindR_analysis_summary.csv")
write_csv(summary_df, summary_file)

# Print final summary
cat("\n===== PathfindR Analysis Summary =====\n")
cat("Projects processed:", nrow(summary_df), "\n")
cat("Projects with enriched terms:", sum(summary_df$Enriched_Terms > 0), "\n")
cat("Total enriched terms:", sum(summary_df$Enriched_Terms), "\n")
cat("Average runtime per project:", round(mean(summary_df$Runtime_Minutes, na.rm = TRUE), 2), "minutes\n")

# Print the summary table
print(summary_df)

cat("\nPathfindR analysis complete for all projects!\n")