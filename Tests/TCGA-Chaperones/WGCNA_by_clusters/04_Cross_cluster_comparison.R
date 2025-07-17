library(readr)
library(dplyr)
library(ggplot2)
library(WGCNA)
library(pheatmap)
library(RColorBrewer)

# Define file paths
results_dir <- "TCGA-Chaperones/WGCNA_by_clusters/results"
chaperone_analysis_dir <- "TCGA-Chaperones/WGCNA_by_clusters/chaperone_analysis"
output_dir <- "TCGA-Chaperones/WGCNA_by_clusters/cross_cluster_comparison"
log_file <- "TCGA-Chaperones/WGCNA_by_clusters/cross_cluster_comparison_log.txt"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
con <- file(log_file, "w")
sink(con, split = TRUE)

# Function to compare modules between clusters within a project
compare_clusters_within_project <- function(project_id) {
  message("\n====== Comparing clusters within ", project_id, " ======")
  
  # Find all cluster directories for this project
  project_results_dir <- file.path(results_dir, project_id)
  cluster_dirs <- list.dirs(project_results_dir, full.names = FALSE, recursive = FALSE)
  cluster_dirs <- cluster_dirs[grepl("^cluster_", cluster_dirs)]
  
  if (length(cluster_dirs) < 2) {
    message("  Less than 2 clusters found for ", project_id, ". Skipping comparison.")
    return(NULL)
  }
  
  # Extract cluster IDs
  cluster_ids <- gsub("cluster_", "", cluster_dirs)
  
  message("  Comparing clusters: ", paste(cluster_ids, collapse = " vs "))
  
  # Load chaperone module data for each cluster
  cluster_data <- list()
  
  for (cluster_id in cluster_ids) {
    chaperone_file <- file.path(chaperone_analysis_dir, project_id, paste0("cluster_", cluster_id), 
                               paste0(project_id, "_cluster_", cluster_id, "_chaperone_modules.csv"))
    
    if (file.exists(chaperone_file)) {
      cluster_data[[cluster_id]] <- read_csv(chaperone_file, show_col_types = FALSE)
    } else {
      message("  Warning: Chaperone file not found for cluster ", cluster_id)
    }
  }
  
  if (length(cluster_data) < 2) {
    message("  Insufficient data for comparison")
    return(NULL)
  }
  
  # Create comparison matrix
  all_genes <- unique(unlist(lapply(cluster_data, function(x) x$base_ensembl)))
  comparison_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(cluster_ids))
  rownames(comparison_matrix) <- all_genes
  colnames(comparison_matrix) <- paste0("Cluster_", cluster_ids)
  
  # Fill the matrix with module assignments
  for (i in seq_along(cluster_data)) {
    cluster_id <- names(cluster_data)[i]
    data <- cluster_data[[cluster_id]]
    
    for (j in seq_len(nrow(data))) {
      gene_id <- data$base_ensembl[j]
      comparison_matrix[gene_id, paste0("Cluster_", cluster_id)] <- data$module_color[j]
    }
  }
  
  # Create output directory for this project
  project_output_dir <- file.path(output_dir, project_id)
  dir.create(project_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save comparison matrix
  comparison_df <- data.frame(
    gene_id = rownames(comparison_matrix),
    comparison_matrix,
    stringsAsFactors = FALSE
  )
  
  write_csv(comparison_df, file.path(project_output_dir, paste0(project_id, "_cluster_comparison.csv")))
  
  # Calculate module consistency
  consistency_stats <- data.frame(
    gene_id = rownames(comparison_matrix),
    n_clusters_found = rowSums(!is.na(comparison_matrix)),
    same_module = apply(comparison_matrix, 1, function(x) {
      non_na <- x[!is.na(x)]
      if (length(non_na) > 1) {
        return(length(unique(non_na)) == 1)
      } else {
        return(NA)
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Add gene symbols
  gene_list <- read_csv("TCGA-Chaperones/gene_list.csv", show_col_types = FALSE)
  consistency_stats <- consistency_stats %>%
    left_join(gene_list, by = c("gene_id" = "ENSEMBL")) %>%
    select(.data$gene_id, .data$Symbol, .data$n_clusters_found, .data$same_module)
  
  write_csv(consistency_stats, file.path(project_output_dir, paste0(project_id, "_gene_consistency.csv")))
  
  # Calculate summary statistics
  genes_in_all_clusters <- sum(consistency_stats$n_clusters_found == length(cluster_ids), na.rm = TRUE)
  genes_consistent <- sum(consistency_stats$same_module, na.rm = TRUE)
  genes_inconsistent <- sum(!consistency_stats$same_module, na.rm = TRUE)
  
  message("  Summary statistics:")
  message("    Genes found in all clusters: ", genes_in_all_clusters)
  message("    Genes with consistent modules: ", genes_consistent)
  message("    Genes with inconsistent modules: ", genes_inconsistent)
  
  # Create visualization
  if (nrow(comparison_df) > 1) {
    # Create heatmap showing module assignments
    # Convert colors to numeric for heatmap
    unique_colors <- unique(c(comparison_matrix))
    unique_colors <- unique_colors[!is.na(unique_colors)]
    
    if (length(unique_colors) > 1) {
      color_to_num <- setNames(seq_along(unique_colors), unique_colors)
      
      numeric_matrix <- apply(comparison_matrix, 2, function(x) {
        sapply(x, function(y) if (is.na(y)) NA else color_to_num[y])
      })
      
      # Create heatmap
      pdf(file.path(project_output_dir, paste0(project_id, "_cluster_comparison_heatmap.pdf")), 
          width = 8, height = 10)
      
      pheatmap(
        numeric_matrix,
        main = paste(project_id, "- Chaperone Module Assignments Across Clusters"),
        color = rainbow(length(unique_colors)),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        fontsize_row = 8,
        fontsize_col = 10,
        na_col = "grey90",
        legend = TRUE
      )
      
      dev.off()
    }
  }
  
  # Return summary
  result <- list(
    project_id = project_id,
    n_clusters = length(cluster_ids),
    cluster_ids = cluster_ids,
    total_genes = nrow(comparison_df),
    genes_in_all_clusters = genes_in_all_clusters,
    genes_consistent = genes_consistent,
    genes_inconsistent = genes_inconsistent,
    consistency_rate = if (genes_consistent + genes_inconsistent > 0) {
      round(genes_consistent / (genes_consistent + genes_inconsistent) * 100, 1)
    } else {
      NA
    }
  )
  
  message("  Comparison completed successfully!")
  message("    Consistency rate: ", result$consistency_rate, "%")
  
  return(result)
}

# Function to compare the same cluster across different projects
compare_same_cluster_across_projects <- function(cluster_id) {
  message("\n====== Comparing Cluster ", cluster_id, " across projects ======")
  
  # Find all projects that have this cluster
  project_dirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
  project_dirs <- project_dirs[project_dirs != ""]
  
  projects_with_cluster <- c()
  cluster_data <- list()
  
  for (project_id in project_dirs) {
    chaperone_file <- file.path(chaperone_analysis_dir, project_id, paste0("cluster_", cluster_id), 
                               paste0(project_id, "_cluster_", cluster_id, "_chaperone_modules.csv"))
    
    if (file.exists(chaperone_file)) {
      projects_with_cluster <- c(projects_with_cluster, project_id)
      cluster_data[[project_id]] <- read_csv(chaperone_file, show_col_types = FALSE)
    }
  }
  
  if (length(projects_with_cluster) < 2) {
    message("  Less than 2 projects have cluster ", cluster_id, ". Skipping comparison.")
    return(NULL)
  }
  
  message("  Comparing cluster ", cluster_id, " across projects: ", paste(projects_with_cluster, collapse = ", "))
  
  # Create comparison matrix
  all_genes <- unique(unlist(lapply(cluster_data, function(x) x$base_ensembl)))
  comparison_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(projects_with_cluster))
  rownames(comparison_matrix) <- all_genes
  colnames(comparison_matrix) <- projects_with_cluster
  
  # Fill the matrix with module assignments
  for (i in seq_along(cluster_data)) {
    project_id <- names(cluster_data)[i]
    data <- cluster_data[[project_id]]
    
    for (j in seq_len(nrow(data))) {
      gene_id <- data$base_ensembl[j]
      comparison_matrix[gene_id, project_id] <- data$module_color[j]
    }
  }
  
  # Create output directory for this cluster
  cluster_output_dir <- file.path(output_dir, paste0("cluster_", cluster_id, "_cross_projects"))
  dir.create(cluster_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save comparison matrix
  comparison_df <- data.frame(
    gene_id = rownames(comparison_matrix),
    comparison_matrix,
    stringsAsFactors = FALSE
  )
  
  write_csv(comparison_df, file.path(cluster_output_dir, paste0("cluster_", cluster_id, "_cross_project_comparison.csv")))
  
  # Calculate module consistency across projects
  consistency_stats <- data.frame(
    gene_id = rownames(comparison_matrix),
    n_projects_found = rowSums(!is.na(comparison_matrix)),
    projects_found = apply(comparison_matrix, 1, function(x) {
      paste(colnames(comparison_matrix)[!is.na(x)], collapse = ", ")
    }),
    stringsAsFactors = FALSE
  )
  
  # Add gene symbols
  gene_list <- read_csv("TCGA-Chaperones/gene_list.csv", show_col_types = FALSE)
  consistency_stats <- consistency_stats %>%
    left_join(gene_list, by = c("gene_id" = "ENSEMBL")) %>%
    select(.data$gene_id, .data$Symbol, .data$n_projects_found, .data$projects_found)
  
  write_csv(consistency_stats, file.path(cluster_output_dir, paste0("cluster_", cluster_id, "_gene_consistency.csv")))
  
  # Calculate summary statistics
  genes_in_all_projects <- sum(consistency_stats$n_projects_found == length(projects_with_cluster), na.rm = TRUE)
  genes_in_most_projects <- sum(consistency_stats$n_projects_found >= ceiling(length(projects_with_cluster) * 0.5), na.rm = TRUE)
  
  message("  Summary statistics:")
  message("    Genes found in all projects: ", genes_in_all_projects)
  message("    Genes found in most projects (≥50%): ", genes_in_most_projects)
  
  # Return summary
  result <- list(
    cluster_id = cluster_id,
    n_projects = length(projects_with_cluster),
    projects = projects_with_cluster,
    total_genes = nrow(comparison_df),
    genes_in_all_projects = genes_in_all_projects,
    genes_in_most_projects = genes_in_most_projects,
    coverage_rate = round(genes_in_all_projects / nrow(comparison_df) * 100, 1)
  )
  
  message("  Cross-project comparison completed successfully!")
  message("    Coverage rate: ", result$coverage_rate, "%")
  
  return(result)
}

# Find all projects with results
project_dirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
project_dirs <- project_dirs[project_dirs != ""]

message("Found results for ", length(project_dirs), " projects")
print(project_dirs)

# Part 1: Compare clusters within each project
message("\n=== PART 1: Within-project cluster comparisons ===")
within_project_results <- list()

for (project_id in project_dirs) {
  result <- try(compare_clusters_within_project(project_id))
  if (!inherits(result, "try-error") && !is.null(result)) {
    within_project_results[[project_id]] <- result
  }
}

# Part 2: Compare same cluster across projects
message("\n=== PART 2: Cross-project cluster comparisons ===")

# Find all unique cluster IDs
all_cluster_ids <- c()
for (project_id in project_dirs) {
  project_results_dir <- file.path(results_dir, project_id)
  cluster_dirs <- list.dirs(project_results_dir, full.names = FALSE, recursive = FALSE)
  cluster_dirs <- cluster_dirs[grepl("^cluster_", cluster_dirs)]
  cluster_ids <- gsub("cluster_", "", cluster_dirs)
  all_cluster_ids <- c(all_cluster_ids, cluster_ids)
}

unique_cluster_ids <- unique(all_cluster_ids)
message("Found unique cluster IDs: ", paste(unique_cluster_ids, collapse = ", "))

cross_project_results <- list()

for (cluster_id in unique_cluster_ids) {
  result <- try(compare_same_cluster_across_projects(cluster_id))
  if (!inherits(result, "try-error") && !is.null(result)) {
    cross_project_results[[cluster_id]] <- result
  }
}

# Generate comprehensive summary
message("\n=== COMPREHENSIVE SUMMARY ===")

# Within-project summary
if (length(within_project_results) > 0) {
  within_project_summary <- data.frame(
    Project = sapply(within_project_results, function(x) x$project_id),
    N_Clusters = sapply(within_project_results, function(x) x$n_clusters),
    Total_Genes = sapply(within_project_results, function(x) x$total_genes),
    Genes_in_All_Clusters = sapply(within_project_results, function(x) x$genes_in_all_clusters),
    Genes_Consistent = sapply(within_project_results, function(x) x$genes_consistent),
    Genes_Inconsistent = sapply(within_project_results, function(x) x$genes_inconsistent),
    Consistency_Rate = sapply(within_project_results, function(x) x$consistency_rate),
    stringsAsFactors = FALSE
  )
  
  write_csv(within_project_summary, file.path(output_dir, "within_project_comparison_summary.csv"))
  
  message("Within-project comparisons:")
  print(within_project_summary)
}

# Cross-project summary
if (length(cross_project_results) > 0) {
  cross_project_summary <- data.frame(
    Cluster_ID = sapply(cross_project_results, function(x) x$cluster_id),
    N_Projects = sapply(cross_project_results, function(x) x$n_projects),
    Total_Genes = sapply(cross_project_results, function(x) x$total_genes),
    Genes_in_All_Projects = sapply(cross_project_results, function(x) x$genes_in_all_projects),
    Genes_in_Most_Projects = sapply(cross_project_results, function(x) x$genes_in_most_projects),
    Coverage_Rate = sapply(cross_project_results, function(x) x$coverage_rate),
    stringsAsFactors = FALSE
  )
  
  write_csv(cross_project_summary, file.path(output_dir, "cross_project_comparison_summary.csv"))
  
  message("\nCross-project comparisons:")
  print(cross_project_summary)
}

# Close log connection
sink()
close(con)

message("\nCross-cluster comparison analysis complete!")
