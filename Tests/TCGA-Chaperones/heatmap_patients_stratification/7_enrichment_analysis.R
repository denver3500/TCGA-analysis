# Load required libraries
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(knitr)

# Set working directory
setwd("/media/windows/BIOINFORMATICS/TCGA/TCGA-analysis/TCGA-Chaperones/heatmap_patients_stratification")

# Define paths
results_dir <- "pathfindR_results"

# Find all pathfindR results files
all_files <- list.files(path = results_dir, 
                        pattern = ".*_pathfindR_results\\.csv$", 
                        recursive = TRUE, 
                        full.names = TRUE)

cat("Found", length(all_files), "pathfindR result files\n")

# Function to extract project name from file path
extract_project <- function(filepath) {
  # Extract project name from path (e.g., TCGA-BLCA from TCGA-BLCA_cluster2_vs_cluster1)
  basename_file <- basename(filepath)
  project_name <- str_extract(basename_file, "TCGA-[A-Z]+")
  comparison <- str_extract(basename_file, "cluster\\d+_vs_cluster\\d+")
  return(data.frame(project = project_name, comparison = comparison))
}

# Process each file and combine results
cat("Processing pathfindR result files...\n")
all_results <- data.frame()

for(file_path in all_files) {
  # Get project info
  project_info <- extract_project(file_path)
  
  # Read the data
  tryCatch({
    data <- read_csv(file_path, show_col_types = FALSE)
    
    # Add project and comparison info
    data$project <- project_info$project
    data$comparison <- project_info$comparison
    
    # Select only essential columns
    data_subset <- data %>%
      select(ID, Term_Description, Fold_Enrichment, lowest_p, project, comparison)
    
    # Add to overall results
    all_results <- bind_rows(all_results, data_subset)
    
  }, error = function(e) {
    cat("Error processing", file_path, ":", e$message, "\n")
  })
}

# Count how many projects each pathway appears in
pathway_summary <- all_results %>%
  group_by(ID, Term_Description) %>%
  summarize(
    times_enriched = n_distinct(project),
    projects = paste(sort(unique(project)), collapse = ", "),
    avg_fold_enrichment = mean(Fold_Enrichment),
    min_p_value = min(lowest_p),
    .groups = "drop"
  ) %>%
  filter(times_enriched > 1) %>%  # Only keep pathways enriched in multiple projects
  arrange(desc(times_enriched), min_p_value)  # Sort by frequency and significance

# Save the results to a CSV file
output_file <- "consolidated_pathfindR_results.csv"
write_csv(pathway_summary, output_file)

# Also create a more detailed file that includes project-specific stats
detailed_summary <- all_results %>%
  filter(ID %in% pathway_summary$ID) %>%  # Only keep pathways enriched in multiple projects
  arrange(ID, project, comparison)

write_csv(detailed_summary, "detailed_pathfindR_results.csv")

# Print a summary table
cat("\nPathways enriched in multiple TCGA projects:\n")
cat("Found", nrow(pathway_summary), "pathways enriched in multiple projects\n")

# Show the top 20 most frequently enriched pathways
top_pathways <- head(pathway_summary, 20)
print(kable(top_pathways, format = "pipe"))

cat("\nResults saved to:", output_file, "and detailed_pathfindR_results.csv\n")

# Create a visualization of the most common pathways
if(require(ggplot2)) {
  # Plot top 20 most frequent pathways
  plot_data <- head(pathway_summary, 20)
  
  p <- ggplot(plot_data, aes(x = reorder(Term_Description, times_enriched), y = times_enriched, fill = -log10(min_p_value))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = "Top Pathways Enriched Across Multiple TCGA Projects",
      x = NULL,
      y = "Number of TCGA Projects",
      fill = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    scale_fill_viridis_c()
  
  # Save the plot
  ggsave("top_enriched_pathways.pdf", p, width = 10, height = 8)
  ggsave("top_enriched_pathways.png", p, width = 10, height = 8, dpi = 300)
  
  cat("Visualizations saved as top_enriched_pathways.pdf and .png\n")
}

# Optional: Create a heatmap showing which pathways are enriched in which projects
if(require(pheatmap) && require(reshape2)) {
  # Get top 30 pathways for heatmap
  top30_pathways <- head(pathway_summary, 30)
  
  # Create a matrix for the heatmap
  heatmap_data <- all_results %>%
    filter(ID %in% top30_pathways$ID) %>%
    mutate(value = -log10(lowest_p)) %>%
    select(ID, Term_Description, project, value)
  
  # Reshape to wide format for heatmap
  heatmap_matrix <- heatmap_data %>%
    distinct(ID, Term_Description, project, value) %>%
    reshape2::dcast(Term_Description ~ project, value.var = "value", fun.aggregate = max, fill = 0)
  
  # Set the rownames to Term_Description and remove that column
  rownames(heatmap_matrix) <- heatmap_matrix$Term_Description
  heatmap_matrix <- heatmap_matrix[, -1]
  
  # Create the heatmap
  pheatmap::pheatmap(
    heatmap_matrix,
    main = "Enriched Pathways Across TCGA Projects",
    filename = "pathways_heatmap.pdf",
    color = colorRampPalette(c("white", "orange", "red", "darkred"))(100),
    display_numbers = TRUE,
    number_format = "%.1f",
    fontsize_row = 8,
    fontsize_col = 10,
    width = 12,
    height = 14
  )
  
  # Also save as PNG
  pheatmap::pheatmap(
    heatmap_matrix,
    main = "Enriched Pathways Across TCGA Projects",
    filename = "pathways_heatmap.png",
    color = colorRampPalette(c("white", "orange", "red", "darkred"))(100),
    display_numbers = TRUE,
    number_format = "%.1f",
    fontsize_row = 8,
    fontsize_col = 10,
    width = 12,
    height = 14
  )
  
  cat("Heatmap visualization saved as pathways_heatmap.pdf and .png\n")
}

cat("\nEnrichment analysis complete!\n")