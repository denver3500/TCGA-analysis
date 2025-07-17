library(WGCNA)
library(readr)
library(dplyr)
library(igraph)

# Suppress R CMD CHECK warnings for variable binding
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("gene_id", "base_ensembl", "Name", "module_color", 
                          "module_number", "is_chaperone", "source", "target", 
                          "correlation", "adjacency", "abs_correlation", 
                          "source_is_chaperone", "target_is_chaperone", 
                          "is_chaperone_edge", "both_chaperones"))
}

# Define file paths
results_dir <- "TCGA-Chaperones/WGCNA_by_clusters_test/results"
output_dir <- "TCGA-Chaperones/WGCNA_by_clusters_test/cytoscape_networks"
log_file <- "TCGA-Chaperones/WGCNA_by_clusters_test/cytoscape_export_log.txt"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
con <- file(log_file, "w")
sink(con, split = TRUE)

# Specify test projects
test_projects <- c("TCGA-BRCA", "TCGA-BLCA")

# Load chaperone gene list
gene_list <- read.csv("TCGA-Chaperones/gene_list.csv", stringsAsFactors = FALSE)
chaperone_genes <- gene_list$ENSEMBL

# Function to export network data for Cytoscape
export_network_for_cytoscape <- function(project_id, cluster_id, correlation_threshold = 0.3) {
  message("\n====== Exporting network for ", project_id, " Cluster ", cluster_id, " ======")
  
  # Load WGCNA workspace
  workspace_file <- file.path(results_dir, project_id, paste0("cluster_", cluster_id), 
                             paste0(project_id, "_cluster_", cluster_id, "_WGCNA_workspace.RData"))
  
  if (!file.exists(workspace_file)) {
    message("  WGCNA workspace not found: ", workspace_file)
    return(NULL)
  }
  
  # Load the workspace
  load(workspace_file)
  message("  Loaded WGCNA workspace")
  
  # Check if required variables exist in workspace
  if (!exists("normalized_data") || !exists("soft_power")) {
    message("  Required variables not found in workspace")
    return(NULL)
  }
  
  # Use the normalized data from the workspace
  wgcna_data <- get("normalized_data")
  soft_power_value <- get("soft_power")
  
  # Load gene modules
  modules_file <- file.path(results_dir, project_id, paste0("cluster_", cluster_id), 
                           paste0(project_id, "_cluster_", cluster_id, "_gene_modules.csv"))
  gene_modules <- read_csv(modules_file, show_col_types = FALSE)
  
  # Extract base ENSEMBL IDs
  extract_base_ensembl <- function(ensembl_id) {
    sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
  }
  
  gene_modules$base_ensembl <- extract_base_ensembl(gene_modules$gene_id)
  
  # Add gene information
  gene_modules <- gene_modules %>%
    left_join(gene_list, by = c("base_ensembl" = "ENSEMBL"))
  
  # Mark chaperone genes
  gene_modules$is_chaperone <- gene_modules$base_ensembl %in% chaperone_genes
  
  # Create output directory for this network
  network_output_dir <- file.path(output_dir, project_id, paste0("cluster_", cluster_id))
  dir.create(network_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Calculate gene correlation matrix
  message("  Calculating gene correlations...")
  gene_cor <- cor(wgcna_data, use = "pairwise.complete.obs")
  
  # Create adjacency matrix using soft power
  message("  Creating adjacency matrix with soft power: ", soft_power_value)
  adjacency <- adjacency(wgcna_data, power = soft_power_value)
  
  # Filter edges by correlation threshold
  message("  Filtering edges with correlation threshold: ", correlation_threshold)
  
  # Get gene pairs above threshold
  edge_list <- data.frame()
  n_genes <- ncol(wgcna_data)
  
  for (i in 1:(n_genes-1)) {
    for (j in (i+1):n_genes) {
      if (abs(gene_cor[i, j]) >= correlation_threshold) {
        gene1 <- colnames(wgcna_data)[i]
        gene2 <- colnames(wgcna_data)[j]
        
        edge_list <- rbind(edge_list, data.frame(
          source = gene1,
          target = gene2,
          correlation = gene_cor[i, j],
          adjacency = adjacency[i, j],
          abs_correlation = abs(gene_cor[i, j])
        ))
      }
    }
  }
  
  message("  Found ", nrow(edge_list), " edges above threshold")
  
  if (nrow(edge_list) == 0) {
    message("  No edges found above threshold. Try lowering the correlation threshold.")
    return(NULL)
  }
  
  # Extract base ENSEMBL for edge list
  edge_list$source_base <- extract_base_ensembl(edge_list$source)
  edge_list$target_base <- extract_base_ensembl(edge_list$target)
  
  # Add gene information to edge list
  edge_list <- edge_list %>%
    left_join(gene_list, by = c("source_base" = "ENSEMBL"), suffix = c("", "_source")) %>%
    left_join(gene_list, by = c("target_base" = "ENSEMBL"), suffix = c("", "_target"))
  
  # Mark chaperone edges
  edge_list$source_is_chaperone <- edge_list$source_base %in% chaperone_genes
  edge_list$target_is_chaperone <- edge_list$target_base %in% chaperone_genes
  edge_list$is_chaperone_edge <- edge_list$source_is_chaperone | edge_list$target_is_chaperone
  edge_list$both_chaperones <- edge_list$source_is_chaperone & edge_list$target_is_chaperone
  
  # Create node list
  all_genes <- unique(c(edge_list$source, edge_list$target))
  node_list <- gene_modules[gene_modules$gene_id %in% all_genes, ]
  
  # Ensure all nodes are included
  missing_nodes <- setdiff(all_genes, node_list$gene_id)
  if (length(missing_nodes) > 0) {
    missing_df <- data.frame(
      gene_id = missing_nodes,
      base_ensembl = extract_base_ensembl(missing_nodes),
      module_number = 0,
      module_color = "grey",
      Name = missing_nodes,
      is_chaperone = extract_base_ensembl(missing_nodes) %in% chaperone_genes
    )
    node_list <- rbind(node_list, missing_df)
  }
  
  # Clean up node list
  node_list <- node_list %>%
    select(gene_id, base_ensembl, Name, module_color, module_number, is_chaperone) %>%
    distinct()
  
  # Add node degree (number of connections)
  node_degrees <- table(c(edge_list$source, edge_list$target))
  node_list$degree <- as.numeric(node_degrees[node_list$gene_id])
  node_list$degree[is.na(node_list$degree)] <- 0
  
  # Export files for Cytoscape
  
  # 1. Edge list (SIF format - Simple Interaction Format)
  sif_file <- file.path(network_output_dir, paste0(project_id, "_cluster_", cluster_id, "_network.sif"))
  sif_data <- data.frame(
    source = edge_list$source,
    interaction = "co-expression",
    target = edge_list$target
  )
  write.table(sif_data, sif_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # 2. Edge attributes
  edge_attr_file <- file.path(network_output_dir, paste0(project_id, "_cluster_", cluster_id, "_edge_attributes.csv"))
  edge_attributes <- edge_list %>%
    select(source, target, correlation, adjacency, abs_correlation, 
           source_is_chaperone, target_is_chaperone, is_chaperone_edge, both_chaperones)
  write_csv(edge_attributes, edge_attr_file)
  
  # 3. Node attributes
  node_attr_file <- file.path(network_output_dir, paste0(project_id, "_cluster_", cluster_id, "_node_attributes.csv"))
  write_csv(node_list, node_attr_file)
  
  # 4. Chaperone-only network
  chaperone_edges <- edge_list[edge_list$is_chaperone_edge, ]
  if (nrow(chaperone_edges) > 0) {
    chaperone_sif_file <- file.path(network_output_dir, paste0(project_id, "_cluster_", cluster_id, "_chaperone_network.sif"))
    chaperone_sif_data <- data.frame(
      source = chaperone_edges$source,
      interaction = "co-expression",
      target = chaperone_edges$target
    )
    write.table(chaperone_sif_data, chaperone_sif_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Chaperone edge attributes
    chaperone_edge_attr_file <- file.path(network_output_dir, paste0(project_id, "_cluster_", cluster_id, "_chaperone_edge_attributes.csv"))
    write_csv(chaperone_edges, chaperone_edge_attr_file)
    
    # Chaperone node attributes
    chaperone_nodes <- unique(c(chaperone_edges$source, chaperone_edges$target))
    chaperone_node_list <- node_list[node_list$gene_id %in% chaperone_nodes, ]
    chaperone_node_attr_file <- file.path(network_output_dir, paste0(project_id, "_cluster_", cluster_id, "_chaperone_node_attributes.csv"))
    write_csv(chaperone_node_list, chaperone_node_attr_file)
  }
  
  # 5. Module-specific networks (for largest modules)
  module_sizes <- table(node_list$module_color)
  large_modules <- names(module_sizes[module_sizes >= 20])  # Modules with at least 20 genes
  
  for (module in large_modules) {
    if (module == "grey") next  # Skip unassigned genes
    
    module_nodes <- node_list$gene_id[node_list$module_color == module]
    module_edges <- edge_list[edge_list$source %in% module_nodes & edge_list$target %in% module_nodes, ]
    
    if (nrow(module_edges) > 0) {
      module_sif_file <- file.path(network_output_dir, paste0(project_id, "_cluster_", cluster_id, "_module_", module, "_network.sif"))
      module_sif_data <- data.frame(
        source = module_edges$source,
        interaction = "co-expression",
        target = module_edges$target
      )
      write.table(module_sif_data, module_sif_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  
  # Create summary
  summary_info <- list(
    project_id = project_id,
    cluster_id = cluster_id,
    total_nodes = nrow(node_list),
    total_edges = nrow(edge_list),
    chaperone_nodes = sum(node_list$is_chaperone),
    chaperone_edges = nrow(chaperone_edges),
    correlation_threshold = correlation_threshold,
    soft_power = soft_power_value,
    n_modules = length(unique(node_list$module_color)),
    output_files = list(
      main_network = basename(sif_file),
      edge_attributes = basename(edge_attr_file),
      node_attributes = basename(node_attr_file),
      chaperone_network = if(nrow(chaperone_edges) > 0) basename(chaperone_sif_file) else NA
    )
  )
  
  message("  Network export completed!")
  message("    Total nodes: ", summary_info$total_nodes)
  message("    Total edges: ", summary_info$total_edges)
  message("    Chaperone nodes: ", summary_info$chaperone_nodes)
  message("    Chaperone edges: ", summary_info$chaperone_edges)
  
  return(summary_info)
}

# Process test projects
message("Exporting networks for Cytoscape analysis...")
message("Test projects: ", paste(test_projects, collapse = ", "))

all_results <- list()

for (project_id in test_projects) {
  message("\n=== Processing project: ", project_id, " ===")
  
  # Find all cluster directories for this project
  project_results_dir <- file.path(results_dir, project_id)
  
  if (!dir.exists(project_results_dir)) {
    message("  No results found for ", project_id)
    next
  }
  
  cluster_dirs <- list.dirs(project_results_dir, full.names = FALSE, recursive = FALSE)
  cluster_dirs <- cluster_dirs[grepl("^cluster_", cluster_dirs)]
  
  if (length(cluster_dirs) == 0) {
    message("  No cluster results found for ", project_id)
    next
  }
  
  # Extract cluster IDs
  cluster_ids <- gsub("cluster_", "", cluster_dirs)
  
  message("  Found clusters: ", paste(cluster_ids, collapse = ", "))
  
  # Process each cluster with different correlation thresholds
  correlation_thresholds <- c(0.3, 0.5, 0.7)
  
  for (cluster_id in cluster_ids) {
    for (threshold in correlation_thresholds) {
      result <- try(export_network_for_cytoscape(project_id, cluster_id, threshold))
      if (!inherits(result, "try-error") && !is.null(result)) {
        result_name <- paste0(project_id, "_cluster_", cluster_id, "_threshold_", threshold)
        all_results[[result_name]] <- result
      }
    }
  }
}

# Generate overall summary
if (length(all_results) > 0) {
  summary_df <- data.frame(
    Project = sapply(all_results, function(x) x$project_id),
    Cluster = sapply(all_results, function(x) x$cluster_id),
    Correlation_Threshold = sapply(all_results, function(x) x$correlation_threshold),
    Total_Nodes = sapply(all_results, function(x) x$total_nodes),
    Total_Edges = sapply(all_results, function(x) x$total_edges),
    Chaperone_Nodes = sapply(all_results, function(x) x$chaperone_nodes),
    Chaperone_Edges = sapply(all_results, function(x) x$chaperone_edges),
    Soft_Power = sapply(all_results, function(x) x$soft_power),
    N_Modules = sapply(all_results, function(x) x$n_modules),
    stringsAsFactors = FALSE
  )
  
  write_csv(summary_df, file.path(output_dir, "cytoscape_export_summary.csv"))
  
  message("\n====== Cytoscape Export Summary ======")
  message("Total networks exported: ", nrow(summary_df))
  print(summary_df)
  
  # Create instruction file
  instructions_file <- file.path(output_dir, "CYTOSCAPE_IMPORT_INSTRUCTIONS.txt")
  
  instructions <- paste0(
    "CYTOSCAPE IMPORT INSTRUCTIONS\n",
    "============================\n\n",
    "Files have been exported for Cytoscape analysis. Here's how to import them:\n\n",
    "1. MAIN NETWORK FILES:\n",
    "   - *.sif files: Network topology (nodes and edges)\n",
    "   - *_edge_attributes.csv: Edge properties (correlation, adjacency, etc.)\n",
    "   - *_node_attributes.csv: Node properties (module, chaperone status, etc.)\n\n",
    "2. SPECIALIZED NETWORKS:\n",
    "   - *_chaperone_network.sif: Networks containing only chaperone gene connections\n",
    "   - *_module_*_network.sif: Networks for specific modules\n\n",
    "3. IMPORT PROCEDURE:\n",
    "   a) Open Cytoscape\n",
    "   b) File > Import > Network from File\n",
    "   c) Select the .sif file you want to analyze\n",
    "   d) File > Import > Table from File\n",
    "   e) Import the corresponding *_node_attributes.csv (as Node Table)\n",
    "   f) Import the corresponding *_edge_attributes.csv (as Edge Table)\n\n",
    "4. CORRELATION THRESHOLDS:\n",
    "   - 0.3: More edges, denser network\n",
    "   - 0.5: Moderate stringency\n",
    "   - 0.7: High stringency, fewer but stronger connections\n\n",
    "5. VISUALIZATION TIPS:\n",
    "   - Color nodes by 'module_color' to see modules\n",
    "   - Size nodes by 'degree' to highlight hubs\n",
    "   - Filter by 'is_chaperone' to focus on chaperone genes\n",
    "   - Edge thickness can represent 'abs_correlation'\n\n",
    "6. RECOMMENDED ANALYSIS:\n",
    "   - Start with threshold 0.5 for balanced network\n",
    "   - Use chaperone-specific networks to focus on your genes of interest\n",
    "   - Compare networks between clusters to identify differences\n",
    "   - Analyze network topology (centrality, clustering coefficient, etc.)\n\n"
  )
  
  writeLines(instructions, instructions_file)
  
} else {
  message("No networks were successfully exported")
}

# Close log connection
sink()
close(con)

message("\nCytoscape network export complete!")
message("Check the 'cytoscape_networks' directory for importable files.")
message("Read CYTOSCAPE_IMPORT_INSTRUCTIONS.txt for detailed import instructions.")
