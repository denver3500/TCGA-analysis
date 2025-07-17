library(WGCNA)
library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# WGCNA parameters
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)

# Define file paths
processed_data_dir <- "TCGA-Chaperones/WGCNA_by_clusters/processed_data"
results_dir <- "TCGA-Chaperones/WGCNA_by_clusters/results"
log_file <- "TCGA-Chaperones/WGCNA_by_clusters/network_construction_log.txt"

# Create results directory
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
con <- file(log_file, "w")
sink(con, split = TRUE)

# Function to construct co-expression network for a cluster
construct_network <- function(project_id, cluster_id) {
  message("\n====== Constructing network for ", project_id, " Cluster ", cluster_id, " ======")
  
  # Load processed data
  data_file <- file.path(processed_data_dir, project_id, paste0(project_id, "_cluster_", cluster_id, "_wgcna_data.rds"))
  
  if (!file.exists(data_file)) {
    message("  Data file not found: ", data_file)
    return(NULL)
  }
  
  wgcna_data <- readRDS(data_file)
  message("  Loaded data: ", nrow(wgcna_data), " samples, ", ncol(wgcna_data), " genes")
  
  # Create output directory for this cluster
  cluster_output_dir <- file.path(results_dir, project_id, paste0("cluster_", cluster_id))
  dir.create(cluster_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Step 1: Check data quality
  message("  Checking data quality...")
  gsg <- goodSamplesGenes(wgcna_data, verbose = 3)
  
  if (!gsg$allOK) {
    message("  Removing outlier genes and samples...")
    if (sum(!gsg$goodGenes) > 0) {
      message("    Removing ", sum(!gsg$goodGenes), " genes")
    }
    if (sum(!gsg$goodSamples) > 0) {
      message("    Removing ", sum(!gsg$goodSamples), " samples")
    }
    wgcna_data <- wgcna_data[gsg$goodSamples, gsg$goodGenes]
  }
  
  # Step 2: Choose soft-thresholding power
  message("  Choosing soft-thresholding power...")
  
  # Test different powers
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  
  # Calculate scale-free topology fit
  sft <- pickSoftThreshold(wgcna_data, powerVector = powers, verbose = 5)
  
  # Plot scale-free topology fit
  pdf(file.path(cluster_output_dir, paste0(project_id, "_cluster_", cluster_id, "_power_selection.pdf")), 
      width = 12, height = 8)
  
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  
  # Scale-free topology fit index
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers, cex = cex1, col = "red")
  abline(h = 0.80, col = "red")
  
  # Mean connectivity
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity", type = "n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  
  dev.off()
  
  # Choose power (first power that gives scale-free R^2 > 0.8, or highest if none do)
  power_threshold <- 0.80
  suitable_powers <- sft$fitIndices[sft$fitIndices[,2] > power_threshold, 1]
  
  if (length(suitable_powers) > 0) {
    soft_power <- min(suitable_powers)
  } else {
    # If no power achieves R^2 > 0.8, choose the one with highest R^2
    soft_power <- sft$fitIndices[which.max(sft$fitIndices[,2]), 1]
  }
  
  message("  Selected soft-thresholding power: ", soft_power)
  message("  Scale-free R^2: ", round(sft$fitIndices[sft$fitIndices[,1] == soft_power, 2], 3))
  
  # Step 3: Construct network and identify modules
  message("  Constructing network and identifying modules...")
  
  # Use automatic module detection
  net <- blockwiseModules(
    wgcna_data,
    power = soft_power,
    TOMType = "unsigned",
    minModuleSize = 30,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamStage = TRUE,
    saveTOMs = TRUE,
    saveTOMFileBase = file.path(cluster_output_dir, paste0(project_id, "_cluster_", cluster_id, "_TOM")),
    verbose = 3
  )
  
  # Convert labels to colors for plotting
  module_colors <- labels2colors(net$colors)
  
  message("  Identified ", length(unique(net$colors)), " modules")
  
  # Save module assignments
  module_assignments <- data.frame(
    gene_id = colnames(wgcna_data),
    module_number = net$colors,
    module_color = module_colors
  )
  
  write_csv(module_assignments, file.path(cluster_output_dir, paste0(project_id, "_cluster_", cluster_id, "_gene_modules.csv")))
  
  # Step 4: Plot dendrograms
  message("  Creating dendrogram plots...")
  
  pdf(file.path(cluster_output_dir, paste0(project_id, "_cluster_", cluster_id, "_dendrograms.pdf")), 
      width = 12, height = 8)
  
  # Plot gene dendrogram and module colors
  plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste(project_id, "- Cluster", cluster_id, "Gene Dendrogram"))
  
  dev.off()
  
  # Step 5: Calculate module eigengenes
  message("  Calculating module eigengenes...")
  
  MEs <- moduleEigengenes(wgcna_data, module_colors)$eigengenes
  
  # Save module eigengenes
  ME_df <- data.frame(
    sample_id = rownames(MEs),
    MEs
  )
  
  write_csv(ME_df, file.path(cluster_output_dir, paste0(project_id, "_cluster_", cluster_id, "_module_eigengenes.csv")))
  
  # Step 6: Create module summary
  module_summary <- data.frame(
    module_color = names(table(module_colors)),
    module_size = as.numeric(table(module_colors))
  ) %>%
    arrange(desc(.data$module_size))
  
  write_csv(module_summary, file.path(cluster_output_dir, paste0(project_id, "_cluster_", cluster_id, "_module_summary.csv")))
  
  # Step 7: Plot module eigengenes heatmap
  message("  Creating module eigengenes heatmap...")
  
  pdf(file.path(cluster_output_dir, paste0(project_id, "_cluster_", cluster_id, "_module_eigengenes_heatmap.pdf")), 
      width = 10, height = 8)
  
  # Calculate dissimilarity of module eigengenes
  ME_dissim <- 1 - cor(MEs)
  
  # Cluster module eigengenes
  ME_tree <- hclust(as.dist(ME_dissim), method = "average")
  
  # Plot the result
  plot(ME_tree, main = paste(project_id, "- Cluster", cluster_id, "Module Eigengenes Clustering"),
       xlab = "", sub = "")
  
  dev.off()
  
  # Step 8: Save workspace
  save(wgcna_data, net, MEs, module_colors, soft_power,
       file = file.path(cluster_output_dir, paste0(project_id, "_cluster_", cluster_id, "_WGCNA_workspace.RData")))
  
  # Return summary
  result <- list(
    project_id = project_id,
    cluster_id = cluster_id,
    n_samples = nrow(wgcna_data),
    n_genes = ncol(wgcna_data),
    soft_power = soft_power,
    scale_free_R2 = sft$fitIndices[sft$fitIndices[,1] == soft_power, 2],
    n_modules = length(unique(net$colors)),
    largest_module_size = max(table(module_colors)),
    output_dir = cluster_output_dir
  )
  
  message("  Network construction completed successfully!")
  message("    Soft power: ", soft_power)
  message("    Scale-free R^2: ", round(result$scale_free_R2, 3))
  message("    Number of modules: ", result$n_modules)
  message("    Largest module size: ", result$largest_module_size)
  
  return(result)
}

# Find all processed data files
processed_projects <- list.dirs(processed_data_dir, full.names = FALSE, recursive = FALSE)
processed_projects <- processed_projects[processed_projects != ""]

message("Found processed data for ", length(processed_projects), " projects")
print(processed_projects)

# Process each project and cluster
all_results <- list()

for (project_id in processed_projects) {
  message("\n=== Processing project: ", project_id, " ===")
  
  # Find all cluster data files for this project
  project_data_dir <- file.path(processed_data_dir, project_id)
  cluster_files <- list.files(project_data_dir, pattern = "*_wgcna_data.rds", full.names = FALSE)
  
  if (length(cluster_files) == 0) {
    message("  No cluster data files found for ", project_id)
    next
  }
  
  # Extract cluster IDs from filenames
  cluster_ids <- gsub(paste0(project_id, "_cluster_"), "", cluster_files)
  cluster_ids <- gsub("_wgcna_data.rds", "", cluster_ids)
  
  message("  Found clusters: ", paste(cluster_ids, collapse = ", "))
  
  # Process each cluster
  for (cluster_id in cluster_ids) {
    result <- try(construct_network(project_id, cluster_id))
    if (!inherits(result, "try-error") && !is.null(result)) {
      all_results[[paste0(project_id, "_cluster_", cluster_id)]] <- result
    }
  }
}

# Generate overall summary
if (length(all_results) > 0) {
  overall_summary <- data.frame(
    Project = sapply(all_results, function(x) x$project_id),
    Cluster = sapply(all_results, function(x) x$cluster_id),
    N_Samples = sapply(all_results, function(x) x$n_samples),
    N_Genes = sapply(all_results, function(x) x$n_genes),
    Soft_Power = sapply(all_results, function(x) x$soft_power),
    Scale_Free_R2 = sapply(all_results, function(x) x$scale_free_R2),
    N_Modules = sapply(all_results, function(x) x$n_modules),
    Largest_Module_Size = sapply(all_results, function(x) x$largest_module_size),
    stringsAsFactors = FALSE
  )
  
  write_csv(overall_summary, file.path(results_dir, "network_construction_summary.csv"))
  
  message("\n====== Overall Network Construction Summary ======")
  message("Total networks constructed: ", nrow(overall_summary))
  message("Average soft power: ", round(mean(overall_summary$Soft_Power), 1))
  message("Average scale-free R^2: ", round(mean(overall_summary$Scale_Free_R2), 3))
  message("Average number of modules: ", round(mean(overall_summary$N_Modules), 1))
  
  print(overall_summary)
} else {
  message("No networks were successfully constructed")
}

# Close log connection
sink()
close(con)

message("\nWGCNA network construction complete!")
