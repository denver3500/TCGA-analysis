# Network generation from WGCNA results
library(WGCNA)
library(readr)
library(dplyr)
library(igraph)
library(RColorBrewer)

# Define directories
wgcna_results_dir <- "TCGA-Chaperones/WGCNA/results"
network_output_dir <- "TCGA-Chaperones/WGCNA/networks"
dir.create(network_output_dir, showWarnings = FALSE, recursive = TRUE)

# Get list of available projects with WGCNA results
project_dirs <- list.dirs(wgcna_results_dir, full.names = FALSE, recursive = FALSE)

# Print available projects
cat("Found", length(project_dirs), "projects with WGCNA results:\n")
print(project_dirs)

# Function to generate network for a project
generate_network <- function(project_id) {
  cat("\n===== Generating network for", project_id, "=====\n")
  
  # Create project output directory
  proj_dir <- file.path(network_output_dir, project_id)
  dir.create(proj_dir, showWarnings = FALSE)
  
  # File paths
  module_file <- file.path(wgcna_results_dir, project_id, paste0(project_id, "_gene_modules.csv"))
  tom_prefix <- file.path(wgcna_results_dir, project_id, paste0(project_id, "_TOM"))
  
  # Load gene module assignments
  if (!file.exists(module_file)) {
    cat("Module assignment file not found for", project_id, "\n")
    return(NULL)
  }
  
  gene_modules <- read_csv(module_file, show_col_types = FALSE)
  
  # Get module colors
  module_colors <- unique(gene_modules$module)
  module_sizes <- table(gene_modules$module)
  
  # Print module information
  cat("Found", length(module_colors), "modules\n")
  print(module_sizes)
  
  # Ask user to select modules for network generation
  cat("\nWhich modules would you like to include in the network?\n")
  cat("Enter module names separated by commas (e.g., 'blue,brown,turquoise')\n")
  cat("For top 3 modules by size, enter 'top3'\n")
  cat("For all modules, enter 'all'\n")
  
  module_selection <- readline(prompt = "Module selection: ")
  
  if (module_selection == "top3") {
    # Get top 3 modules by size
    top_modules <- names(sort(module_sizes, decreasing = TRUE))[1:min(3, length(module_colors))]
    cat("Selected top 3 modules:", paste(top_modules, collapse = ", "), "\n")
  } else if (module_selection == "all") {
    # Use all modules
    top_modules <- module_colors
    cat("Selected all modules\n")
  } else {
    # Parse user selection
    top_modules <- strsplit(module_selection, ",")[[1]]
    top_modules <- trimws(top_modules)
    # Check if all selected modules exist
    invalid_modules <- top_modules[!top_modules %in% module_colors]
    if (length(invalid_modules) > 0) {
      cat("Warning: The following modules do not exist:", paste(invalid_modules, collapse = ", "), "\n")
      top_modules <- top_modules[top_modules %in% module_colors]
    }
    cat("Selected modules:", paste(top_modules, collapse = ", "), "\n")
  }
  
  if (length(top_modules) == 0) {
    cat("No valid modules selected. Exiting.\n")
    return(NULL)
  }
  
  # Load TOM data
  # Find all TOM files
  tom_files <- list.files(dirname(tom_prefix), pattern = paste0(basename(tom_prefix), "-block.[0-9]+.RData"), full.names = TRUE)
  
  if (length(tom_files) == 0) {
    cat("No TOM files found for", project_id, "\n")
    return(NULL)
  }
  
  cat("Found", length(tom_files), "TOM files\n")
  
  # Create empty list to store edge lists from all blocks
  all_edges <- list()
  
  # Process each TOM block
  for (tom_file in tom_files) {
    cat("Processing", basename(tom_file), "\n")
    
    # Load the TOM for this block
    load(tom_file)
    
    # Check if TOM exists in the loaded data
    if (!exists("TOM")) {
      cat("Warning: TOM object not found in", basename(tom_file), "\n")
      next
    }
    
    # Get block number from filename
    block_num <- as.numeric(gsub(".*-block\\.([0-9]+)\\.RData", "\\1", tom_file))
    cat("Block", block_num, "has", nrow(TOM), "genes\n")
    
    # Check gene ID matching between module assignment and TOM matrix
    tom_genes <- rownames(TOM)
    module_genes_all <- gene_modules$gene_id
    
    # Print sample gene IDs from both sources
    cat("First few rownames in TOM:", paste(head(tom_genes), collapse=", "), "\n")
    cat("First few gene IDs in module file:", paste(head(module_genes_all), collapse=", "), "\n")
    
    # Check overlap between TOM and module gene IDs
    overlap_count <- sum(module_genes_all %in% tom_genes)
    cat("Number of genes that overlap between module file and TOM:", overlap_count, 
        "(", round(overlap_count/length(module_genes_all)*100, 1), "%)\n")
    
    # Create a gene ID mapping function
    map_gene_ids <- function(module_genes) {
      if (overlap_count >= length(module_genes_all) * 0.1) {
        # Regular matching by name if sufficient overlap
        return(match(module_genes, tom_genes))
      } else {
        cat("WARNING: Very few gene IDs match between module file and TOM matrix\n")
        
        # Try matching by removing version numbers from Ensembl IDs
        module_genes_base <- gsub("\\.[0-9]+$", "", module_genes)
        tom_genes_base <- gsub("\\.[0-9]+$", "", tom_genes)
        base_indices <- match(module_genes_base, tom_genes_base)
        
        base_overlap_count <- sum(!is.na(base_indices))
        cat("After removing version numbers, found", base_overlap_count, "matches\n")
        
        if (base_overlap_count >= length(module_genes) * 0.1) {
          return(base_indices)
        }
        
        # If still poor match, try positional matching if dimensions are compatible
        if (length(tom_genes) == length(module_genes_all)) {
          cat("Using positional matching as a last resort\n")
          # Create a mapping between positions in the two lists
          id_mapping <- data.frame(
            module_index = 1:length(module_genes_all),
            tom_index = 1:length(tom_genes)
          )
          
          # Get indices of selected module genes in the mapping
          module_indices <- match(module_genes, module_genes_all)
          
          # Map to TOM indices
          return(id_mapping$tom_index[module_indices])
        } else {
          cat("WARNING: Cannot create a reliable mapping between gene IDs\n")
          return(NULL)
        }
      }
    }
    
    # For each selected module, get genes
    for (mod in top_modules) {
      # Get genes in this module
      module_genes <- gene_modules$gene_id[gene_modules$module == mod]
      cat("Module", mod, "has", length(module_genes), "genes\n")
      
      # Map gene IDs between module file and TOM matrix
      gene_indices <- map_gene_ids(module_genes)
      
      # Remove NA indices
      gene_indices <- gene_indices[!is.na(gene_indices)]
      
      cat("Found", length(gene_indices), "genes from module", mod, "in block", block_num, "\n")
      
      if (length(gene_indices) < 2) {
        # Not enough genes from this module in this block
        cat("Not enough genes from module", mod, "in block", block_num, "to create edges\n")
        next
      }
      
      # Extract submatrix for this module
      mod_tom <- TOM[gene_indices, gene_indices]
      
      # Convert to edge list
      if (length(gene_indices) > 1) {  # Need at least 2 genes to create edges
        gene_names <- rownames(TOM)[gene_indices]
        
        # Create edge list
        edges <- data.frame()
        for (i in 1:(length(gene_indices)-1)) {
          for (j in (i+1):length(gene_indices)) {
            weight <- mod_tom[i, j]
            
            # Only include edges with weight above threshold
            # Lower the threshold to include more connections
            if (weight > 0.001) {  # Lowered threshold from 0.01
              edge <- data.frame(
                source = gene_names[i],
                target = gene_names[j],
                weight = weight,
                module1 = mod,
                module2 = mod
              )
              edges <- rbind(edges, edge)
            }
          }
        }
        
        if (nrow(edges) > 0) {
          cat("Generated", nrow(edges), "edges for module", mod, "in block", block_num, "\n")
          all_edges[[length(all_edges) + 1]] <- edges
        } else {
          cat("No edges above threshold for module", mod, "in block", block_num, "\n")
        }
      }
    }
  }
  
  # Combine all edges
  if (length(all_edges) > 0) {
    edge_list <- do.call(rbind, all_edges)
    
    # Sort by weight
    edge_list <- edge_list[order(-edge_list$weight),]
    
    # Save edge list
    edge_file <- file.path(proj_dir, paste0(project_id, "_edge_list.tsv"))
    write_tsv(edge_list, edge_file)
    
    cat("Saved edge list with", nrow(edge_list), "edges to", edge_file, "\n")
    
    # Identify hub genes (genes with highest connectivity)
    source_counts <- table(edge_list$source)
    target_counts <- table(edge_list$target)
    gene_connections <- table(c(names(source_counts), names(target_counts)))
    
    # Get top hub genes
    top_hub_genes <- sort(gene_connections, decreasing = TRUE)[1:min(20, length(gene_connections))]
    
    # Save hub genes
    hub_genes <- data.frame(
      gene = names(top_hub_genes),
      connections = as.numeric(top_hub_genes)
    )
    
    hub_file <- file.path(proj_dir, paste0(project_id, "_hub_genes.csv"))
    write_csv(hub_genes, hub_file)
    
    cat("Saved top", nrow(hub_genes), "hub genes to", hub_file, "\n")
    
    # Create simplified network for visualization (limit edges if too many)
    if (nrow(edge_list) > 1000) {
      cat("Edge list is large (", nrow(edge_list), "edges). Creating simplified network for visualization.\n")
      
      # Option 1: Take top N edges by weight
      simplified_edges <- edge_list[1:1000,]
      
      simplified_file <- file.path(proj_dir, paste0(project_id, "_simplified_network.tsv"))
      write_tsv(simplified_edges, simplified_file)
      cat("Saved simplified network with", nrow(simplified_edges), "edges to", simplified_file, "\n")
    }
    
    return(list(
      project_id = project_id,
      edge_count = nrow(edge_list),
      modules = top_modules,
      hub_genes = hub_genes$gene[1:min(10, nrow(hub_genes))]
    ))
  } else {
    cat("No edges found for the selected modules\n")
    
    # Add alternative method if no edges found with TOM approach
    cat("Trying alternative method using direct gene correlations...\n")
    
    # Try loading expression data
    expr_file <- file.path("TCGA-Chaperones/WGCNA/processed_data", paste0(project_id, "_high_expr_wgcna.rds"))
    if (file.exists(expr_file)) {
      expr_data <- readRDS(expr_file)
      cat("Loaded expression data:", nrow(expr_data), "samples x", ncol(expr_data), "genes\n")
      
      # Create edges for each module using correlation
      alt_edges <- data.frame()
      
      for (mod in top_modules) {
        cat("Processing module", mod, "with direct correlation\n")
        module_genes <- gene_modules$gene_id[gene_modules$module == mod]
        
        # Find these genes in expression data
        gene_cols <- match(module_genes, colnames(expr_data))
        gene_cols <- gene_cols[!is.na(gene_cols)]
        
        if (length(gene_cols) > 1) {
          # Calculate correlation
          expr_subset <- expr_data[, gene_cols]
          cor_mat <- cor(expr_subset)
          
          # Create edges for highly correlated genes
          for (i in 1:(length(gene_cols)-1)) {
            for (j in (i+1):length(gene_cols)) {
              if (abs(cor_mat[i,j]) > 0.6) { # Only strong correlations
                edge <- data.frame(
                  source = colnames(expr_subset)[i],
                  target = colnames(expr_subset)[j],
                  weight = abs(cor_mat[i,j]),
                  module1 = mod,
                  module2 = mod
                )
                alt_edges <- rbind(alt_edges, edge)
              }
            }
          }
        }
      }
      
      if (nrow(alt_edges) > 0) {
        cat("Generated", nrow(alt_edges), "edges using correlation method\n")
        
        # Save alternative edge list
        alt_edge_file <- file.path(proj_dir, paste0(project_id, "_correlation_edges.tsv"))
        write_tsv(alt_edges, alt_edge_file)
        
        # Generate hub genes
        source_counts <- table(alt_edges$source)
        target_counts <- table(alt_edges$target)
        gene_connections <- table(c(names(source_counts), names(target_counts)))
        top_hub_genes <- sort(gene_connections, decreasing = TRUE)[1:min(20, length(gene_connections))]
        
        hub_genes <- data.frame(
          gene = names(top_hub_genes),
          connections = as.numeric(top_hub_genes)
        )
        
        alt_hub_file <- file.path(proj_dir, paste0(project_id, "_correlation_hub_genes.csv"))
        write_csv(hub_genes, alt_hub_file)
        
        return(list(
          project_id = project_id,
          edge_count = nrow(alt_edges),
          modules = top_modules,
          hub_genes = hub_genes$gene[1:min(10, nrow(hub_genes))],
          method = "correlation"
        ))
      }
    }
    
    return(NULL)
  }
}

# Let user select a project
selected_project_id <- readline(prompt = "Enter a project ID to generate network (e.g., TCGA-THYM): ")

# Check if the selected project exists
if (selected_project_id %in% project_dirs) {
  # Generate network for the selected project
  result <- generate_network(selected_project_id)
  if (!is.null(result)) {
    cat("\nNetwork generation complete for", selected_project_id, "\n")
    cat("Results saved to", file.path(network_output_dir, selected_project_id), "\n")
    
    cat("\nTop hub genes:\n")
    print(result$hub_genes)
    
    method_note <- if(is.null(result$method)) "" else " (generated using correlation method)"
    
    cat("\nTo visualize this network in Cytoscape:\n")
    cat("1. Open Cytoscape\n")
    cat("2. Import > Network > File...\n")
    cat("3. Select the edge list file", method_note, "\n")
    cat("4. Use 'weight' column for edge width/opacity\n")
    cat("5. Use 'module1' column for node colors\n")
  }
} else {
  cat("Project", selected_project_id, "not found. Available projects:", paste(project_dirs, collapse=", "), "\n")
}