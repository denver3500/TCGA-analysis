library(readr)
library(dplyr)
library(ggplot2)
library(WGCNA)
library(pheatmap)
library(RColorBrewer)

# Load the gene list
gene_list <- read.csv("TCGA-Chaperones/gene_list.csv", stringsAsFactors = FALSE)

# Extract ENSEMBL IDs from gene list
chaperone_genes <- gene_list$ENSEMBL

# Define file paths
results_dir <- "TCGA-Chaperones/WGCNA_by_clusters/results"
output_dir <- "TCGA-Chaperones/WGCNA_by_clusters/chaperone_analysis"
log_file <- "TCGA-Chaperones/WGCNA_by_clusters/chaperone_analysis_log.txt"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
con <- file(log_file, "w")
sink(con, split = TRUE)

# Function to analyze chaperone genes in modules
analyze_chaperone_modules <- function(project_id, cluster_id) {
  message("\n====== Analyzing chaperone genes for ", project_id, " Cluster ", cluster_id, " ======")
  
  # File paths
  cluster_results_dir <- file.path(results_dir, project_id, paste0("cluster_", cluster_id))
  gene_modules_file <- file.path(cluster_results_dir, paste0(project_id, "_cluster_", cluster_id, "_gene_modules.csv"))
  
  if (!file.exists(gene_modules_file)) {
    message("  Gene modules file not found: ", gene_modules_file)
    return(NULL)
  }
  
  # Load gene modules
  gene_modules <- read_csv(gene_modules_file, show_col_types = FALSE)
  
  # Extract base ENSEMBL IDs (remove version numbers)
  extract_base_ensembl <- function(ensembl_id) {
    sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
  }
  
  gene_modules$base_ensembl <- extract_base_ensembl(gene_modules$gene_id)
  
  # Find chaperone genes in modules
  chaperone_modules <- gene_modules[gene_modules$base_ensembl %in% chaperone_genes, ]
  
  if (nrow(chaperone_modules) == 0) {
    message("  No chaperone genes found in modules")
    return(NULL)
  }
  
  message("  Found ", nrow(chaperone_modules), " chaperone genes in modules")
  
  # Add gene names from gene list
  chaperone_modules <- chaperone_modules %>%
    left_join(gene_list, by = c("base_ensembl" = "ENSEMBL")) %>%
    select(.data$gene_id, .data$base_ensembl, .data$Symbol, .data$module_number, .data$module_color)
  
  # Count chaperone genes per module
  module_counts <- chaperone_modules %>%
    group_by(.data$module_color) %>%
    summarise(
      n_chaperones = n(),
      chaperone_genes = paste(.data$Symbol, collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(desc(.data$n_chaperones))
  
  message("  Chaperone genes distribution across modules:")
  print(module_counts)
  
  # Calculate percentages
  total_chaperones <- nrow(chaperone_modules)
  module_counts$percentage <- round(module_counts$n_chaperones / total_chaperones * 100, 1)
  
  # Identify dominant module
  dominant_module <- module_counts$module_color[1]
  dominant_count <- module_counts$n_chaperones[1]
  dominant_percentage <- module_counts$percentage[1]
  
  message("  Dominant module: ", dominant_module, " (", dominant_count, " genes, ", dominant_percentage, "%)")
  
  # Create output directory for this analysis
  analysis_output_dir <- file.path(output_dir, project_id, paste0("cluster_", cluster_id))
  dir.create(analysis_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save detailed results
  write_csv(chaperone_modules, file.path(analysis_output_dir, paste0(project_id, "_cluster_", cluster_id, "_chaperone_modules.csv")))
  write_csv(module_counts, file.path(analysis_output_dir, paste0(project_id, "_cluster_", cluster_id, "_module_counts.csv")))
  
  # Create visualization
  if (nrow(module_counts) > 1) {
    # Bar plot of chaperone distribution
    p1 <- ggplot(module_counts, aes(x = reorder(.data$module_color, .data$n_chaperones), y = .data$n_chaperones)) +
      geom_bar(stat = "identity", fill = module_counts$module_color) +
      coord_flip() +
      labs(title = paste(project_id, "- Cluster", cluster_id, "Chaperone Gene Distribution"),
           x = "Module Color",
           y = "Number of Chaperone Genes") +
      theme_minimal() +
      theme(axis.text.y = element_text(color = module_counts$module_color[order(module_counts$n_chaperones)]))
    
    ggsave(file.path(analysis_output_dir, paste0(project_id, "_cluster_", cluster_id, "_chaperone_distribution.pdf")), 
           p1, width = 10, height = 6)
    
    # Pie chart
    p2 <- ggplot(module_counts, aes(x = "", y = .data$n_chaperones, fill = .data$module_color)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      scale_fill_identity() +
      labs(title = paste(project_id, "- Cluster", cluster_id, "Chaperone Module Distribution"),
           fill = "Module") +
      theme_void() +
      theme(legend.position = "right")
    
    ggsave(file.path(analysis_output_dir, paste0(project_id, "_cluster_", cluster_id, "_chaperone_pie.pdf")), 
           p2, width = 8, height = 6)
  }
  
  # Return summary
  result <- list(
    project_id = project_id,
    cluster_id = cluster_id,
    total_chaperones = total_chaperones,
    total_modules = nrow(module_counts),
    dominant_module = dominant_module,
    dominant_count = dominant_count,
    dominant_percentage = dominant_percentage,
    module_distribution = module_counts,
    chaperone_details = chaperone_modules
  )
  
  message("  Analysis completed successfully!")
  
  return(result)
}

# Find all cluster result directories
project_dirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
project_dirs <- project_dirs[project_dirs != ""]

message("Found results for ", length(project_dirs), " projects")
print(project_dirs)

# Process each project and cluster
all_results <- list()

for (project_id in project_dirs) {
  message("\n=== Processing project: ", project_id, " ===")
  
  # Find all cluster directories for this project
  project_results_dir <- file.path(results_dir, project_id)
  cluster_dirs <- list.dirs(project_results_dir, full.names = FALSE, recursive = FALSE)
  cluster_dirs <- cluster_dirs[grepl("^cluster_", cluster_dirs)]
  
  if (length(cluster_dirs) == 0) {
    message("  No cluster results found for ", project_id)
    next
  }
  
  # Extract cluster IDs
  cluster_ids <- gsub("cluster_", "", cluster_dirs)
  
  message("  Found clusters: ", paste(cluster_ids, collapse = ", "))
  
  # Process each cluster
  for (cluster_id in cluster_ids) {
    result <- try(analyze_chaperone_modules(project_id, cluster_id))
    if (!inherits(result, "try-error") && !is.null(result)) {
      all_results[[paste0(project_id, "_cluster_", cluster_id)]] <- result
    }
  }
}

# Generate comprehensive summary
if (length(all_results) > 0) {
  message("\n====== Comprehensive Summary ======")
  
  # Summary statistics
  summary_stats <- data.frame(
    Project = sapply(all_results, function(x) x$project_id),
    Cluster = sapply(all_results, function(x) x$cluster_id),
    Total_Chaperones = sapply(all_results, function(x) x$total_chaperones),
    Total_Modules = sapply(all_results, function(x) x$total_modules),
    Dominant_Module = sapply(all_results, function(x) x$dominant_module),
    Dominant_Count = sapply(all_results, function(x) x$dominant_count),
    Dominant_Percentage = sapply(all_results, function(x) x$dominant_percentage),
    stringsAsFactors = FALSE
  )
  
  write_csv(summary_stats, file.path(output_dir, "chaperone_analysis_summary.csv"))
  
  # Create cross-project comparison
  message("Cross-project chaperone module analysis:")
  print(summary_stats)
  
  # Analyze clustering consistency
  message("\nClustering consistency analysis:")
  clustering_consistency <- summary_stats %>%
    group_by(.data$Project) %>%
    summarise(
      n_clusters = n(),
      avg_chaperones_per_cluster = mean(.data$Total_Chaperones),
      avg_dominant_percentage = mean(.data$Dominant_Percentage),
      consistent_clustering = length(unique(.data$Dominant_Module)) == 1,
      .groups = "drop"
    )
  
  print(clustering_consistency)
  write_csv(clustering_consistency, file.path(output_dir, "clustering_consistency.csv"))
  
  # Create heatmap of dominant module percentages
  heatmap_data <- summary_stats %>%
    select(.data$Project, .data$Cluster, .data$Dominant_Percentage) %>%
    tidyr::pivot_wider(names_from = .data$Cluster, values_from = .data$Dominant_Percentage, values_fill = 0) %>%
    as.data.frame()
  
  rownames(heatmap_data) <- heatmap_data$Project
  heatmap_data <- heatmap_data[, -1]
  
  if (ncol(heatmap_data) > 1 && nrow(heatmap_data) > 1) {
    pdf(file.path(output_dir, "chaperone_clustering_heatmap.pdf"), width = 10, height = 8)
    pheatmap(
      heatmap_data,
      main = "Chaperone Gene Clustering Consistency\n(% in Dominant Module)",
      color = colorRampPalette(c("white", "red"))(100),
      display_numbers = TRUE,
      fontsize_number = 10,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      angle_col = 45
    )
    dev.off()
  }
  
  # Overall statistics
  message("\n====== Overall Statistics ======")
  message("Total analyses: ", nrow(summary_stats))
  message("Projects analyzed: ", length(unique(summary_stats$Project)))
  message("Average chaperones per cluster: ", round(mean(summary_stats$Total_Chaperones), 1))
  message("Average dominant module percentage: ", round(mean(summary_stats$Dominant_Percentage), 1), "%")
  message("Clusters with >50% chaperones in dominant module: ", sum(summary_stats$Dominant_Percentage > 50))
  message("Clusters with >75% chaperones in dominant module: ", sum(summary_stats$Dominant_Percentage > 75))
  
} else {
  message("No analyses were successfully completed")
}

# Close log connection
sink()
close(con)

message("\nChaperone module analysis complete!")
