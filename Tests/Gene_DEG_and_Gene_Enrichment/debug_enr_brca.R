library(pathfindR)
library(readr)
library(stringr)
library(dplyr)

# Define the specific project and gene to analyze
project_id <- "TCGA-BRCA"
gene_name <- "BRCA1"  # Updated for your CSV file

# Define just the two combinations to test
gene_sets_options <- c("KEGG", "GO-All")
pin_name_options <- c("KEGG")  # Only using KEGG PIN

# Set up logging to both console and file
log_file <- paste0("TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/pathfindR_combinations_", gene_name, "_", project_id, ".txt")
con <- file(log_file, "w")
sink(con, split=TRUE)  # Write output to both console and file

message("====== PATHFINDR PARAMETER COMBINATIONS ANALYSIS START ======")
message(paste("Project:", project_id))
message(paste("Gene:", gene_name))  
message(paste("Testing", length(gene_sets_options) * length(pin_name_options), "parameter combinations"))

# Define paths for your specific CSV file
input_file <- "TCGA-BRCA_BRCA1_DEG_results_debug.csv"
results_dir <- paste0("TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/", gene_name, "_pathfindR_combinations")

# Create results directory
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Check if the CSV file exists
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# Read the CSV file
message(paste("Reading input file:", input_file))
deg_data <- read_csv(input_file, show_col_types = FALSE)

# First, check if we need to add gene names
if ("gene_name" %in% colnames(deg_data)) {
  message("Gene names already present in CSV file")
} else if ("gene_ensembl" %in% colnames(deg_data)) {
  message("Need to add gene names - looking up from RDS file")
  
  # Find and load the corresponding RDS file
  rds_dir <- "TCGA-Chaperones/rds"
  rds_pattern <- paste0("^", project_id, "_.*\\.rds$")
  rds_files <- list.files(rds_dir, pattern = rds_pattern, full.names = TRUE)
  
  if (length(rds_files) == 0) {
    stop("Could not find RDS file for project ", project_id, " in directory ", rds_dir)
  }
  
  # Extract base ENSEMBL ID without version function
  extract_base_ensembl <- function(ensembl_id) {
    # Split by period and take the first part (base ID without version)
    return(sapply(strsplit(ensembl_id, "\\."), function(x) x[1]))
  }
  
  # Load the RDS file to extract gene metadata
  se_object <- readRDS(rds_files[1])
  
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
    
    # Extract base ENSEMBL IDs from the DEG data
    deg_data$ensembl_base <- extract_base_ensembl(deg_data$gene_ensembl)
    
    # Use base R merge instead of dplyr
    mapping_subset <- ensembl_to_gene_map[, c("ensembl_base", "gene_name")]
    deg_data <- merge(deg_data, mapping_subset, by="ensembl_base", all.x=TRUE)
    
    # Report how many IDs were successfully mapped
    mapped_count <- sum(!is.na(deg_data$gene_name))
    message("Successfully mapped ", mapped_count, " out of ", nrow(deg_data), " entries")
  } else {
    stop("RDS file does not have expected gene metadata structure")
  }
}

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

# Loop through just the two combinations
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
message("\nTop parameter combinations by number of enriched terms:")
print(sorted_summary[, c("gene_sets", "pin_name", "num_terms", "runtime_seconds")])

# Close the log connection
sink() 
close(con)

message("Analysis complete! Results saved to:", results_dir)
message("Log file saved to:", log_file)
message("Summary file saved to:", file.path(results_dir, "combinations_summary.csv"))