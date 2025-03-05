# Script to create a CSV summary from the saved RDS file

# Load the saved RDS file
rds_file <- file.path("TCGA-WGCNA", "TCGA_gene_expression_analysis_results.rds")
cat("Loading data from:", rds_file, "\n")
combined_results <- readRDS(rds_file)

# Extract match results
match_results <- combined_results$match_results
cat("Extracted match results for", length(match_results), "projects\n")

# Find projects with matches
projects_with_matches <- names(match_results)[sapply(match_results, function(x) x$has_matches)]
cat("Found", length(projects_with_matches), "projects with gene matches\n")

if (length(projects_with_matches) > 0) {
  # Create a data frame to store the information
  match_summary_df <- data.frame(
    Project = character(),
    Match_Count = integer(),
    Matched_Genes = character(),
    stringsAsFactors = FALSE
  )
  
  # Fill the data frame with match information
  for (project in projects_with_matches) {
    match_info <- match_results[[project]]
    
    # Add row to data frame
    match_summary_df <- rbind(match_summary_df, data.frame(
      Project = project,
      Match_Count = match_info$match_count,
      Matched_Genes = paste(match_info$matches, collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by number of matches (descending)
  match_summary_df <- match_summary_df[order(match_summary_df$Match_Count, decreasing = TRUE), ]
  
  # Write to CSV
  csv_output_file <- file.path("analysis", "project_gene_variance_matches.csv")
  write.csv(match_summary_df, file = csv_output_file, row.names = FALSE)
  cat("Match summary saved to:", csv_output_file, "\n")
  
  # Print top 5 projects with most matches
  cat("\nTop 5 projects with most matches:\n")
  print(head(match_summary_df, 5))
} else {
  cat("\nNo projects with matches found.\n")
}