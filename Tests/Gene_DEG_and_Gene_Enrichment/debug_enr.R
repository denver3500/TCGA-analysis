library(pathfindR)
library(readr)
library(stringr)
library(dplyr)

# Define the specific project and gene to analyze
project_id <- "TCGA-BRCA"
gene_name <- "HSP90AA1"

# Define all parameter combinations to test
gene_sets_options <- c("KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC", "GO-MF", "cell_markers")
pin_name_options <- c("Biogrid", "STRING", "GeneMania", "IntAct", "KEGG")

# Set up logging to both console and file
log_file <- paste0("TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/pathfindR_combinations_", gene_name, "_", project_id, ".txt")
con <- file(log_file, "w")
sink(con, split=TRUE)  # Write output to both console and file

message("====== PATHFINDR PARAMETER COMBINATIONS ANALYSIS START ======")
message(paste("Project:", project_id))
message(paste("Gene:", gene_name))
message(paste("Testing", length(gene_sets_options) * length(pin_name_options), "parameter combinations"))

# Define paths
project_path <- file.path("TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/DEG_results", project_id)
deg_file <- file.path(project_path, paste0(gene_name, "_DEG_results.csv"))
results_dir <- file.path(project_path, paste0(gene_name, "_pathfindR_combinations"))

# Create results directory
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Check if DEG results file exists
if (!file.exists(deg_file)) {
  stop(paste("DEG results file not found:", deg_file))
}

# Read the DEG file
deg_data <- read_csv(deg_file, show_col_types = FALSE)

# Prepare input dataframe for pathfindR
if ("gene_name" %in% colnames(deg_data) && "logFC" %in% colnames(deg_data)) {
  if ("FDR" %in% colnames(deg_data)) {
    input_df <- data.frame(
      Gene_symbol = deg_data$gene_name,
      logFC = deg_data$logFC,
      FDR_p = deg_data$FDR
    )
  } else if ("PValue" %in% colnames(deg_data)) {
    input_df <- data.frame(
      Gene_symbol = deg_data$gene_name,
      logFC = deg_data$logFC,
      FDR_p = deg_data$PValue
    )
  } else {
    stop("Neither FDR nor PValue columns found in DEG file")
  }
} else {
  stop("Required columns not found in DEG file")
}

# Remove rows with NA in Gene_symbol
input_df <- input_df[!is.na(input_df$Gene_symbol), ]
message(paste("Input data has", nrow(input_df), "genes with valid symbols"))

# Create a summary table to track results
summary_table <- data.frame(
  gene_sets = character(),
  pin_name = character(),
  success = logical(),
  num_terms = numeric(),
  error = character(),
  runtime_seconds = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all combinations
combination_count <- 0
total_combinations <- length(gene_sets_options) * length(pin_name_options)

for (gene_sets in gene_sets_options) {
  for (pin_name in pin_name_options) {
    combination_count <- combination_count + 1
    
    # Create a unique name for this combination
    combo_name <- paste0(gene_sets, "_", pin_name)
    message("\n============================================")
    message(paste("Running combination", combination_count, "of", total_combinations, ":", combo_name))
    
    # Create directory for this combination
    combo_dir <- file.path(results_dir, combo_name)
    if (dir.exists(combo_dir)) {
      unlink(combo_dir, recursive = TRUE)
    }
    dir.create(combo_dir, showWarnings = FALSE)
    
    # Initialize result variables
    success <- FALSE
    num_terms <- 0
    error_message <- NA
    runtime <- NA
    
    # Run pathfindR with the current parameter combination
    tryCatch({
      message(paste("Starting pathfindR run with gene_sets =", gene_sets, "and pin_name =", pin_name))
      start_time <- Sys.time()
      
      pathfinder_results <- run_pathfindR(
        input_df,
        output_dir = combo_dir,
        gene_sets = gene_sets,
        pin_name = pin_name,
        min_gset_size = 5,    # Using more permissive settings
        max_gset_size = 500,
        p_val_threshold = 0.05,
        enrichment_threshold = 0,
        silent_option = TRUE
      )
      
      end_time <- Sys.time()
      runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      # Check if we have results
      if (is.data.frame(pathfinder_results) && nrow(pathfinder_results) > 0) {
        success <- TRUE
        num_terms <- nrow(pathfinder_results)
        
        # Save results
        result_file <- file.path(combo_dir, paste0(gene_name, "_", combo_name, "_results.csv"))
        write_csv(pathfinder_results, result_file)
        
        message(paste("Success! Found", num_terms, "enriched terms in", round(runtime, 2), "seconds"))
        message("Top 5 terms:")
        if (nrow(pathfinder_results) > 0) {
          print(head(pathfinder_results[, c("ID", "Term_Description", "Adjusted_P_value")], 5))
        }
      } else {
        message("No enriched terms found")
        success <- TRUE  # Still consider it a success, just with 0 terms
        num_terms <- 0
      }
    }, error = function(e) {
      message(paste("Error:", e$message))
      error_message <- e$message
      success <- FALSE
    })
    
    # Add to summary table
    summary_table <- rbind(summary_table, data.frame(
      gene_sets = gene_sets,
      pin_name = pin_name,
      success = success,
      num_terms = num_terms,
      error = ifelse(is.na(error_message), "None", error_message),
      runtime_seconds = ifelse(is.na(runtime), NA, round(runtime, 2)),
      stringsAsFactors = FALSE
    ))
    
    # Save the updated summary table after each combination
    write_csv(summary_table, file.path(results_dir, "combinations_summary.csv"))
    
    message(paste("Completed combination", combination_count, "of", total_combinations))
  }
}

# Print overall summary
message("\n====== OVERALL SUMMARY ======")
message(paste("Total combinations tested:", nrow(summary_table)))
message(paste("Successful runs:", sum(summary_table$success)))
message(paste("Runs with errors:", sum(!summary_table$success)))
message(paste("Combinations with enriched terms:", sum(summary_table$num_terms > 0)))

# Sort summary by number of terms found
sorted_summary <- summary_table[order(-summary_table$num_terms), ]
message("\nTop 5 parameter combinations by number of enriched terms:")
print(head(sorted_summary[, c("gene_sets", "pin_name", "num_terms", "runtime_seconds")], 5))

# Close the log connection
sink() 
close(con)

message("Analysis complete! Results saved to:", results_dir)
message("Log file saved to:", log_file)
message("Summary file saved to:", file.path(results_dir, "combinations_summary.csv"))