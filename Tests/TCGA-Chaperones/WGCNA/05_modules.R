# Load the gene list
gene_list <- read.csv("TCGA-Chaperones/gene_list.csv", stringsAsFactors = FALSE)

# Extract ENSEMBL IDs from gene list
chaperone_genes <- gene_list$Name

# Get all project folders
results_dir <- "TCGA-Chaperones/WGCNA/results"
project_folders <- list.dirs(results_dir, recursive = FALSE, full.names = FALSE)
project_folders <- project_folders[grepl("^TCGA-", project_folders)]

# Initialize a list to store results
all_results <- list()

# Loop through each project
for(project in project_folders) {
  cat("\n", strrep("=", 60), "\n")
  cat("Analyzing project:", project, "\n")
  cat(strrep("=", 60), "\n")
  
  # Load the gene modules file for this project
  modules_file <- file.path(results_dir, project, paste0(project, "_gene_modules.csv"))
  
  if(file.exists(modules_file)) {
    gene_modules <- read.csv(modules_file, stringsAsFactors = FALSE)
    
    # Check which chaperone genes are in the modules file and their assignments
    chaperone_modules <- gene_modules[gene_modules$gene_id %in% chaperone_genes, ]
    
    # Display the results
    print("Chaperone genes and their module assignments:")
    print(chaperone_modules)
    
    # Count genes per module
    module_counts <- table(chaperone_modules$module)
    print("\nNumber of chaperone genes per module:")
    print(module_counts)
    
    # Calculate percentages
    total_found <- nrow(chaperone_modules)
    print(paste("\nTotal chaperone genes found in modules:", total_found, "out of", length(chaperone_genes)))
    
    if(total_found > 0) {
      percentages <- round(module_counts / total_found * 100, 1)
      print("\nPercentage of chaperone genes per module:")
      print(percentages)
      
      # Identify the dominant module
      dominant_module <- names(which.max(module_counts))
      dominant_percentage <- max(percentages)
      print(paste("\nDominant module:", dominant_module, "with", dominant_percentage, "% of chaperone genes"))
      
      # Store results for summary
      all_results[[project]] <- list(
        total_found = total_found,
        module_counts = module_counts,
        percentages = percentages,
        dominant_module = dominant_module,
        dominant_percentage = dominant_percentage,
        chaperone_modules = chaperone_modules
      )
    } else {
      all_results[[project]] <- list(
        total_found = 0,
        module_counts = NULL,
        percentages = NULL,
        dominant_module = NA,
        dominant_percentage = 0,
        chaperone_modules = data.frame()
      )
    }
  } else {
    cat("Gene modules file not found for", project, "\n")
  }
}

# Summary across all projects
cat("\n", strrep("=", 60), "\n")
cat("SUMMARY ACROSS ALL PROJECTS\n")
cat(strrep("=", 60), "\n")

summary_df <- data.frame(
  Project = names(all_results),
  Total_Found = sapply(all_results, function(x) x$total_found),
  Dominant_Module = sapply(all_results, function(x) ifelse(is.na(x$dominant_module), "None", x$dominant_module)),
  Dominant_Percentage = sapply(all_results, function(x) x$dominant_percentage),
  stringsAsFactors = FALSE
)

print(summary_df)

# Create a comprehensive module assignment table
cat("\nModule assignments across all projects:\n")
all_assignments <- data.frame()

for(project in names(all_results)) {
  if(nrow(all_results[[project]]$chaperone_modules) > 0) {
    temp_df <- all_results[[project]]$chaperone_modules
    temp_df$Project <- project
    all_assignments <- rbind(all_assignments, temp_df[, c("Project", "gene_id", "module")])
  }
}

if(nrow(all_assignments) > 0) {
  # Reshape to wide format for easy comparison
  library(tidyr)
  wide_assignments <- all_assignments %>%
    pivot_wider(names_from = Project, values_from = module, values_fill = "Not_Found")
  
  print(wide_assignments)
  
  # Save results
  write.csv(summary_df, "TCGA-Chaperones/WGCNA/chaperone_module_summary.csv", row.names = FALSE)
  write.csv(all_assignments, "TCGA-Chaperones/WGCNA/chaperone_module_assignments.csv", row.names = FALSE)
  write.csv(wide_assignments, "TCGA-Chaperones/WGCNA/chaperone_module_assignments_wide.csv", row.names = FALSE)
}