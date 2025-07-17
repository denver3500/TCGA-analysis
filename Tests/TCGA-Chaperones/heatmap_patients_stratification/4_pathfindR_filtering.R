# Load required libraries
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# Set working directory
setwd("/media/windows/BIOINFORMATICS/TCGA/TCGA-analysis/TCGA-Chaperones/heatmap_patients_stratification")

# Read the gene list
gene_list <- read_csv("../gene_list.csv")
chaperone_genes <- gene_list$Name

# Function to check if any chaperone genes are present in the Up_regulated or Down_regulated columns
contains_chaperone_genes <- function(up_genes, down_genes, target_genes) {
  # Split gene lists by comma and trim whitespace
  up_list <- if(!is.na(up_genes)) str_split(up_genes, ",")[[1]] %>% str_trim() else character(0)
  down_list <- if(!is.na(down_genes)) str_split(down_genes, ",")[[1]] %>% str_trim() else character(0)
  
  # Combine all genes
  all_genes <- c(up_list, down_list)
  
  # Check if any target genes are present
  any(target_genes %in% all_genes)
}

# Function to process a single pathfindR results file
process_pathfindR_file <- function(file_path, output_dir, project_name) {
  # Read the CSV file
  data <- read_csv(file_path, show_col_types = FALSE)
  
  # Filter rows that contain chaperone genes
  filtered_data <- data %>%
    rowwise() %>%
    filter(contains_chaperone_genes(Up_regulated, Down_regulated, chaperone_genes)) %>%
    ungroup()
  
  # Only create output file if there are filtered results
  if (nrow(filtered_data) > 0) {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create output filename
    output_file <- file.path(output_dir, paste0(project_name, "_chaperone_filtered.csv"))
    
    # Write filtered data
    write_csv(filtered_data, output_file)
    
    cat("Processed", project_name, ":", nrow(filtered_data), "pathways containing chaperone genes\n")
    return(nrow(filtered_data))
  } else {
    cat("No chaperone genes found in", project_name, "\n")
    return(0)
  }
}

# Main processing function
process_all_pathfindR_files <- function() {
  # Set paths
  pathfindR_dir <- "pathfindR_results"
  output_dir <- "pathfindR_results/chaperone_filtered"
  
  # Get all subdirectories (TCGA projects)
  project_dirs <- list.dirs(pathfindR_dir, full.names = TRUE, recursive = FALSE)
  project_dirs <- project_dirs[!grepl("chaperone_filtered", project_dirs)]  # Exclude output directory
  
  # Initialize summary
  summary_results <- data.frame(
    Project = character(),
    Pathways_with_chaperones = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each project directory
  for (project_dir in project_dirs) {
    project_name <- basename(project_dir)
    
    # Find CSV files in the directory
    csv_files <- list.files(project_dir, pattern = "\\.csv$", full.names = TRUE)
    
    if (length(csv_files) > 0) {
      # Process the first CSV file found (assuming one per directory)
      csv_file <- csv_files[1]
      
      # Process the file
      pathway_count <- process_pathfindR_file(csv_file, output_dir, project_name)
      
      # Add to summary
      summary_results <- rbind(summary_results, 
                               data.frame(Project = project_name, 
                                         Pathways_with_chaperones = pathway_count))
    } else {
      cat("No CSV files found in", project_name, "\n")
    }
  }
  
  # Write summary
  write_csv(summary_results, file.path(output_dir, "filtering_summary.csv"))
  
  # Print summary
  cat("\n=== SUMMARY ===\n")
  print(summary_results)
  cat("\nTotal projects processed:", nrow(summary_results), "\n")
  cat("Projects with chaperone pathways:", sum(summary_results$Pathways_with_chaperones > 0), "\n")
  
  return(summary_results)
}

# Run the analysis
cat("Starting pathfindR filtering for chaperone genes...\n")
cat("Chaperone genes to search for:", paste(chaperone_genes, collapse = ", "), "\n\n")

results_summary <- process_all_pathfindR_files()

cat("\nFiltering complete! Check the 'pathfindR_results/chaperone_filtered' directory for results.\n")