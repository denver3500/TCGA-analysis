library(TCGAbiolinks)

# Load the saved RDS file
rds_file <- file.path("TCGA-WGCNA", "TCGA_gene_expression_analysis_results.rds")
cat("Loading data from:", rds_file, "\n")
combined_results <- readRDS(rds_file)

# Get TCGA project information 
projects_info <- TCGAbiolinks:::getGDCprojects()
cat("Loaded information for", nrow(projects_info), "projects\n")

# Extract match results
match_results <- combined_results$match_results
projects_with_matches <- names(match_results)[sapply(match_results, function(x) x$has_matches)]

if (length(projects_with_matches) > 0) {
  # Create our data frame
  match_summary_df <- data.frame(
    project_id = gsub("_transcriptomic_exp\\.rds$", "", projects_with_matches),
    Match_Count = sapply(match_results[projects_with_matches], function(x) x$match_count),
    Matched_Genes = sapply(match_results[projects_with_matches], function(x) paste(x$matches, collapse = ";")),
    stringsAsFactors = FALSE
  )
  
  # Add project names by direct lookup in projects_info
  match_summary_df$project_name <- sapply(match_summary_df$project_id, function(pid) {
    # Look for the exact ID
    matches <- projects_info$name[projects_info$id == pid]
    if(length(matches) > 0) {
      return(matches[1])
    }
    
    # Try with TCGA- prefix if not found and if not already prefixed
    if (!grepl("^TCGA-", pid)) {
      tcga_pid <- paste0("TCGA-", pid)
      matches <- projects_info$name[projects_info$id == tcga_pid]
      if(length(matches) > 0) {
        return(matches[1])
      }
    }
    
    # Return unknown if no match found
    return(paste0("Unknown (", pid, ")"))
  })
  
  # Reorder columns to desired format
  final_df <- match_summary_df[, c("project_name", "Match_Count", "Matched_Genes", "project_id")]
  
  # Create output directory
  if (!dir.exists("analysis")) {
    dir.create("analysis")
  }
  
  # Write to CSV
  csv_output_file <- file.path("analysis", "project_gene_variance_matches.csv")
  write.csv(final_df, file = csv_output_file, row.names = FALSE)
  
} else {
  cat("No projects with matches found.\n")
}