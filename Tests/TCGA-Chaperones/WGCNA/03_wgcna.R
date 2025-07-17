# Consolidated WGCNA analysis and network generation
# Following https://bioinformaticsworkbook.org/tutorials/wgcna.html

# Load required libraries
library(WGCNA)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)
library(igraph)
library(RColorBrewer)

# Enable multi-threading for WGCNA
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 4) # Adjust based on your system

# Define directories
wgcna_data_dir <- "TCGA-Chaperones/WGCNA/processed_data"
annotation_dir <- "TCGA-Chaperones/WGCNA/gene_annotation"
output_dir <- "TCGA-Chaperones/WGCNA/results"
network_dir <- "TCGA-Chaperones/WGCNA/networks"

# Create directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(network_dir, showWarnings = FALSE, recursive = TRUE)

# Get list of available projects
high_expr_files <- list.files(wgcna_data_dir, pattern = "_high_expr_wgcna\\.rds$", full.names = TRUE)
project_ids <- gsub(".*/(TCGA-[A-Z]+)_high_expr_wgcna\\.rds$", "\\1", high_expr_files)

# Print available projects
cat("Found", length(project_ids), "projects for WGCNA analysis:\n")
print(project_ids)

# Function to annotate expression data and handle duplicates
annotate_expression_data <- function(expr_data, project_id, proj_dir) {
  # Load annotation file
  annotation_file <- file.path(annotation_dir, paste0(project_id, "_wgcna_gene_mapping.csv"))
  
  # Create mapping to store original and annotated IDs
  gene_mapping <- data.frame(
    original_id = colnames(expr_data),
    final_id = colnames(expr_data),
    stringsAsFactors = FALSE
  )
  
  if (!file.exists(annotation_file)) {
    cat("Warning: No annotation file found for", project_id, "\n")
    cat("Using original gene IDs\n")
  } else {
    gene_data <- read_csv(annotation_file, show_col_types = FALSE)
    
    # Match gene IDs from expression data to annotation
    genes_in_data <- colnames(expr_data)
    matched_indices <- match(genes_in_data, gene_data$ensembl_id)
    
    # If some genes don't match, try matching with base Ensembl IDs (without version)
    if (any(is.na(matched_indices))) {
      base_genes_in_data <- gsub("\\.[0-9]+$", "", genes_in_data)
      for (i in which(is.na(matched_indices))) {
        potential_match <- which(gene_data$ensembl_base_id == base_genes_in_data[i])
        if (length(potential_match) > 0) {
          matched_indices[i] <- potential_match[1]
        }
      }
    }
    
    # Create a mapping from original IDs to gene symbols
    gene_symbols <- gene_data$gene_symbol[matched_indices]
    
    # Replace NA with original IDs
    gene_symbols[is.na(gene_symbols)] <- genes_in_data[is.na(gene_symbols)]
    
    # Update the mapping dataframe
    gene_mapping$final_id <- gene_symbols
    
    # Handle duplicate gene symbols by keeping the one with higher average expression
    if (any(duplicated(gene_symbols))) {
      cat("Found", sum(duplicated(gene_symbols)), "duplicate gene symbols\n")
      cat("Keeping genes with higher average expression\n")
      
      # Add mean expression to the mapping
      gene_mapping$mean_expr <- colMeans(expr_data, na.rm = TRUE)
      
      # Find duplicates and keep highest expression for each symbol
      gene_mapping <- gene_mapping %>%
        arrange(final_id, desc(mean_expr)) %>%
        distinct(final_id, .keep_all = TRUE)
      
      # Subset expression data to keep only non-duplicate genes
      expr_data <- expr_data[, match(gene_mapping$original_id, colnames(expr_data))]
      
      # Update column names to gene symbols
      colnames(expr_data) <- gene_mapping$final_id
    } else {
      # Just rename columns if no duplicates
      colnames(expr_data) <- gene_symbols
    }
  }
  
  # Save the gene ID mapping for later reference
  write_csv(gene_mapping, file.path(proj_dir, paste0(project_id, "_gene_id_mapping.csv")))
  cat("Saved gene ID mapping to", file.path(proj_dir, paste0(project_id, "_gene_id_mapping.csv")), "\n")
  
  return(expr_data)
}

# Combined function for WGCNA analysis and network generation
run_wgcna_workflow <- function(project_id, generate_network = TRUE) {
  cat("\n===== Running WGCNA for", project_id, "=====\n")
  
  # Create project output directories
  proj_dir <- file.path(output_dir, project_id)
  network_proj_dir <- file.path(network_dir, project_id)
  dir.create(proj_dir, showWarnings = FALSE)
  dir.create(network_proj_dir, showWarnings = FALSE)
  
  # Load expression data
  expr_file <- file.path(wgcna_data_dir, paste0(project_id, "_high_expr_wgcna.rds"))
  expr_data <- readRDS(expr_file)
  
  # Check data dimensions
  cat("Loaded expression data:", nrow(expr_data), "samples x", ncol(expr_data), "genes\n")
  
  # Save original column names before annotation (for reference in network generation)
  original_gene_ids <- colnames(expr_data)
  
  # Annotate expression data and get gene ID mapping
  expr_data <- annotate_expression_data(expr_data, project_id, proj_dir)
  cat("After annotation and duplicate removal:", nrow(expr_data), "samples x", ncol(expr_data), "genes\n")
  
  # Check for missing values
  if (sum(is.na(expr_data)) > 0) {
    cat("Warning: Found", sum(is.na(expr_data)), "missing values. Replacing with column means.\n")
    expr_data <- impute::impute.knn(expr_data)$data
  }
  
  # Filter genes with low variance
  gene_vars <- apply(expr_data, 2, var)
  high_var_genes <- gene_vars > quantile(gene_vars, 0.25)  # Keep top 75% variable genes
  expr_data <- expr_data[, high_var_genes]
  cat("After filtering low-variance genes:", nrow(expr_data), "samples x", ncol(expr_data), "genes\n")
  
  # Step 1: Choose soft-thresholding power
  # Choose power range based on dataset size
  powers <- c(1:10, seq(from = 12, to = 30, by = 2))
  
  # Calculate scale-free topology fit indices for various powers
  cat("Analyzing network topology for various soft-thresholding powers...\n")
  sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 0)
  
  # Create diagnostic plots
  par(mfrow = c(1, 2))
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       main = paste(project_id, "- Scale-Free Topology Fit"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = 0.7, col = "red")
  abline(h = 0.8, col = "red", lty = "dashed")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       main = paste(project_id, "- Mean Connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels = powers, cex = 0.7, col = "red")
  
  # Save the diagnostic plots
  dev.copy(pdf, file.path(proj_dir, paste0(project_id, "_power_selection.pdf")), width = 10, height = 6)
  dev.off()
  
  # Reset plot layout
  par(mfrow = c(1, 1))
  
  # Get user input for power selection
  cat("\nBased on the scale-free topology fit, choose a soft-thresholding power.\n")
  cat("Look for the lowest power where the scale-free topology fit index curve flattens out.\n")
  cat("Recommended value from analysis:", sft$powerEstimate, "\n")
  cat("Available powers:", paste(powers, collapse = ", "), "\n")
  
  # Get user input
  picked_power <- as.numeric(readline(prompt = "Enter your picked soft-thresholding power: "))
  
  # Check if the input is valid
  while (!picked_power %in% powers) {
    cat("Invalid power value. Please choose from:", paste(powers, collapse = ", "), "\n")
    picked_power <- as.numeric(readline(prompt = "Enter your picked soft-thresholding power: "))
  }
  
  cat("Using soft-thresholding power:", picked_power, "\n")
  
  # Step 2: Build the network and identify modules
  cat("Constructing network and detecting modules...\n")
  
  net <- blockwiseModules(expr_data, 
                         power = picked_power,
                         TOMType = "unsigned", 
                         minModuleSize = 30,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = file.path(proj_dir, paste0(project_id, "_TOM")),
                         verbose = 3)
  
  # Step 3: Convert numeric labels to color labels
  module_colors <- labels2colors(net$colors)
  
  # Get module sizes
  module_sizes <- table(module_colors)
  cat("\nDetected", length(module_sizes), "modules\n")
  print(module_sizes)
  
  # Get the number of blocks
  block_count <- length(net$dendrograms)
  cat("\nNetwork was split into", block_count, "blocks\n")
  
  # Module visualization
  pdf(file.path(proj_dir, paste0(project_id, "_module_heatmap.pdf")), width = 10, height = 10)
  # Prepare data for color-coded module display
  modules <- as.data.frame(table(module_colors))
  colnames(modules) <- c("Module", "GeneCount")
  modules <- modules[order(-modules$GeneCount),]
  
  # Plot module sizes as a bar chart
  barplot(modules$GeneCount, 
          names.arg = modules$Module,
          col = as.character(modules$Module),
          main = paste(project_id, "- Module Sizes"),
          ylab = "Number of genes",
          las = 2)
  
  # Display the top modules
  head(modules, 10)
  dev.off()
  
  # Create separate dendrogram plots for each block if there are multiple blocks
  if (block_count > 0) {
    pdf(file.path(proj_dir, paste0(project_id, "_block_dendrograms.pdf")), width = 12, height = 8)
    
    # If there are many blocks, limit the number displayed
    blocks_to_plot <- min(block_count, 6)
    
    # Set up the plotting layout
    layout_rows <- ceiling(sqrt(blocks_to_plot))
    layout_cols <- ceiling(blocks_to_plot / layout_rows)
    par(mfrow = c(layout_rows, layout_cols))
    
    for (i in 1:blocks_to_plot) {
      # Check the structure of blockGenes to ensure proper access
      cat("Processing block", i, "\n")
      
      # Correctly access genes in this block - blockGenes is a list, not a vector
      if (is.list(net$blockGenes)) {
        # For newer WGCNA versions where blockGenes is a list
        if (length(net$blockGenes) >= i) {
          block_genes <- net$blockGenes[[i]]
        } else {
          cat("Block", i, "not found in blockGenes list\n")
          next
        }
      } else {
        # For older versions that might use a vector with block assignments
        tryCatch({
          block_genes <- which(as.vector(net$blockGenes) == i)
        }, error = function(e) {
          cat("Error accessing block genes:", conditionMessage(e), "\n")
          cat("Skipping block", i, "\n")
          return(NULL)
        })
      }
      
      # Skip if no genes found in this block
      if (length(block_genes) == 0) {
        cat("No genes found in block", i, "\n")
        next
      }
      
      cat("Found", length(block_genes), "genes in block", i, "\n")
      
      # Get the module colors for these genes
      block_colors <- module_colors[block_genes]
      
      # Plot the block dendrogram if it exists
      if (i <= length(net$dendrograms)) {
        plotDendroAndColors(net$dendrograms[[i]], colors = block_colors,
                          groupLabels = paste("Block", i),
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = paste("Block", i, "-", length(block_genes), "genes"))
      }
    }
    dev.off()
  }
  
  # Step 4: Relate modules to each other
  # Calculate eigengenes
  ME_list <- moduleEigengenes(expr_data, colors = module_colors)
  MEs <- ME_list$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  ME_diss <- 1 - cor(MEs)
  
  # Cluster module eigengenes
  ME_tree <- hclust(as.dist(ME_diss), method = "average")
  
  # Plot the module relationship dendrogram
  pdf(file.path(proj_dir, paste0(project_id, "_module_relationships.pdf")), width = 10, height = 7)
  plot(ME_tree, main = paste(project_id, "- Module Relationships"),
       xlab = "", sub = "", cex = 0.7)
  dev.off()
  
  # Step 5: Save module membership
  # Create data frame with gene module assignments
  gene_module_df <- data.frame(
    gene_id = colnames(expr_data),
    module = module_colors
  )
  
  # Save module assignments
  write_csv(gene_module_df, file.path(proj_dir, paste0(project_id, "_gene_modules.csv")))
  
  # Save module eigengenes
  write_csv(as.data.frame(MEs), file.path(proj_dir, paste0(project_id, "_module_eigengenes.csv")))
  
  # Step 6: Generate network if requested
  if (generate_network) {
    cat("\n===== Generating network for", project_id, "=====\n")
    
    # Ask user to select modules for network generation
    cat("\nWhich modules would you like to include in the network?\n")
    cat("Enter module names separated by commas (e.g., 'blue,brown,turquoise')\n")
    cat("For top 3 modules by size, enter 'top3'\n")
    cat("For all modules, enter 'all'\n")
    
    module_selection <- readline(prompt = "Module selection: ")
    
    if (module_selection == "top3") {
      # Get top 3 modules by size
      top_modules <- names(sort(module_sizes, decreasing = TRUE))[1:min(3, length(unique(module_colors)))]
      cat("Selected top 3 modules:", paste(top_modules, collapse = ", "), "\n")
    } else if (module_selection == "all") {
      # Use all modules
      top_modules <- unique(module_colors)
      cat("Selected all modules\n")
    } else {
      # Parse user selection
      top_modules <- strsplit(module_selection, ",")[[1]]
      top_modules <- trimws(top_modules)
      # Check if all selected modules exist
      invalid_modules <- top_modules[!top_modules %in% unique(module_colors)]
      if (length(invalid_modules) > 0) {
        cat("Warning: The following modules do not exist:", paste(invalid_modules, collapse = ", "), "\n")
        top_modules <- top_modules[top_modules %in% unique(module_colors)]
      }
      cat("Selected modules:", paste(top_modules, collapse = ", "), "\n")
    }
    
    if (length(top_modules) == 0) {
      cat("No valid modules selected. Skipping network generation.\n")
      return(list(
        project_id = project_id,
        module_colors = module_colors,
        module_sizes = module_sizes,
        eigengenes = MEs,
        power = picked_power,
        network_generated = FALSE
      ))
    }
    
    # We have two approaches for generating networks:
    # 1. Using TOM files (saved during blockwiseModules)
    # 2. Using direct correlations between genes
    
    # Load gene ID mapping to ensure we're using consistent IDs
    id_mapping <- read_csv(file.path(proj_dir, paste0(project_id, "_gene_id_mapping.csv")), 
                           show_col_types = FALSE)
    
    # Create empty list to store edge lists from all blocks/methods
    all_edges <- list()
    
    # Approach 1: Try using TOM files first
    tom_files <- list.files(proj_dir, pattern = paste0(project_id, "_TOM-block.[0-9]+.RData"), full.names = TRUE)
    
    if (length(tom_files) > 0) {
      cat("Found", length(tom_files), "TOM files\n")
      
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
        
        # For each selected module, get genes
        for (mod in top_modules) {
          # Get genes in this module
          module_gene_ids <- gene_module_df$gene_id[gene_module_df$module == mod]
          cat("Module", mod, "has", length(module_gene_ids), "genes\n")
          
          # Find these genes in the TOM
          # We need to map back to original IDs that were used in the TOM
          
          # Look at TOM rownames - they might be original or annotated IDs
          tom_genes <- rownames(TOM)
          sample_tom <- head(tom_genes, 5)
          cat("Sample TOM rownames:", paste(sample_tom, collapse=", "), "\n")
          
          # Let's check if the rownames in TOM are the original IDs or annotated IDs
          if (all(sample_tom %in% id_mapping$original_id)) {
            # TOM uses original IDs
            cat("TOM uses original gene IDs\n")
            # Map module genes back to original IDs
            original_ids <- id_mapping$original_id[match(module_gene_ids, id_mapping$final_id)]
            gene_indices <- match(original_ids, tom_genes)
          } else if (all(sample_tom %in% id_mapping$final_id)) {
            # TOM uses annotated IDs
            cat("TOM uses annotated gene IDs\n")
            gene_indices <- match(module_gene_ids, tom_genes)
          } else {
            # No match - try using positional indices if ordering is preserved
            cat("WARNING: Cannot match gene IDs between modules and TOM\n")
            cat("Will try a correlation-based approach instead\n")
            next
          }
          
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
                # A lowered threshold to include more connections
                if (weight > 0.001) {  
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
    } else {
      cat("No TOM files found. Using correlation-based approach.\n")
    }
    
    # If no edges were found with TOM, use correlation approach
    if (length(all_edges) == 0) {
      cat("Using correlation-based approach for network generation...\n")
      
      # Process each module using correlations directly from expression data
      for (mod in top_modules) {
        cat("Processing module", mod, "with direct correlation\n")
        
        # Get genes in this module
        module_gene_ids <- gene_module_df$gene_id[gene_module_df$module == mod]
        
        if (length(module_gene_ids) < 2) {
          cat("Module", mod, "has too few genes (", length(module_gene_ids), ") to create edges\n")
          next
        }
        
        # Get expression data for these genes
        genes_in_expr <- match(module_gene_ids, colnames(expr_data))
        genes_in_expr <- genes_in_expr[!is.na(genes_in_expr)]
        
        if (length(genes_in_expr) > 1) {
          # Calculate correlation
          expr_subset <- expr_data[, genes_in_expr]
          cor_mat <- cor(expr_subset)
          
          # Create edges for highly correlated genes
          edges <- data.frame()
          for (i in 1:(ncol(expr_subset)-1)) {
            for (j in (i+1):ncol(expr_subset)) {
              if (abs(cor_mat[i,j]) > 0.6) { # Only strong correlations
                edge <- data.frame(
                  source = colnames(expr_subset)[i],
                  target = colnames(expr_subset)[j],
                  weight = abs(cor_mat[i,j]),
                  module1 = mod,
                  module2 = mod
                )
                edges <- rbind(edges, edge)
              }
            }
          }
          
          if (nrow(edges) > 0) {
            cat("Generated", nrow(edges), "edges for module", mod, "using correlation\n")
            all_edges[[length(all_edges) + 1]] <- edges
          } else {
            cat("No strong correlations found for module", mod, "\n")
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
      edge_file <- file.path(network_proj_dir, paste0(project_id, "_edge_list.tsv"))
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
      
      hub_file <- file.path(network_proj_dir, paste0(project_id, "_hub_genes.csv"))
      write_csv(hub_genes, hub_file)
      
      cat("Saved top", nrow(hub_genes), "hub genes to", hub_file, "\n")
      
      # Create simplified network for visualization (limit edges if too many)
      if (nrow(edge_list) > 1000) {
        cat("Edge list is large (", nrow(edge_list), "edges). Creating simplified network for visualization.\n")
        
        # Take top N edges by weight
        simplified_edges <- edge_list[1:1000,]
        
        simplified_file <- file.path(network_proj_dir, paste0(project_id, "_simplified_network.tsv"))
        write_tsv(simplified_edges, simplified_file)
        cat("Saved simplified network with", nrow(simplified_edges), "edges to", simplified_file, "\n")
      }
      
      # Print visualization instructions
      cat("\nNetwork generation complete for", project_id, "\n")
      cat("Results saved to", network_proj_dir, "\n")
      
      cat("\nTop hub genes:\n")
      print(hub_genes$gene[1:min(10, nrow(hub_genes))])
      
      cat("\nTo visualize this network in Cytoscape:\n")
      cat("1. Open Cytoscape\n")
      cat("2. Import > Network > File...\n")
      cat("3. Select the edge list file\n")
      cat("4. Use 'weight' column for edge width/opacity\n")
      cat("5. Use 'module1' column for node colors\n")
      
      network_generated <- TRUE
    } else {
      cat("No edges were generated for any module\n")
      network_generated <- FALSE
    }
    
    # Return results
    return(list(
      project_id = project_id,
      module_colors = module_colors,
      module_sizes = module_sizes,
      eigengenes = MEs,
      power = picked_power,
      network_generated = network_generated,
      hub_genes = if(exists("hub_genes")) hub_genes$gene else NULL
    ))
  } else {
    # Return just WGCNA results if network generation was not requested
    return(list(
      project_id = project_id,
      module_colors = module_colors,
      module_sizes = module_sizes,
      eigengenes = MEs,
      power = picked_power,
      network_generated = FALSE
    ))
  }
}

# Ask user what they want to do
cat("\nWhat would you like to do?\n")
cat("1: Run WGCNA analysis only\n")
cat("2: Run network generation only (requires previous WGCNA results)\n")
cat("3: Run full workflow (WGCNA analysis + network generation)\n")
choice <- as.numeric(readline(prompt = "Enter your choice (1-3): "))

if (choice == 1) {
  # WGCNA analysis only
  selected_project_id <- readline(prompt = "Enter a project ID to analyze (e.g., TCGA-THYM): ")
  if (selected_project_id %in% project_ids) {
    result <- run_wgcna_workflow(selected_project_id, generate_network = FALSE)
    cat("\nWGCNA analysis complete for", selected_project_id, "\n")
    cat("Results saved to", file.path(output_dir, selected_project_id), "\n")
  } else {
    cat("Project", selected_project_id, "not found. Available projects:", paste(project_ids, collapse=", "), "\n")
  }
} else if (choice == 2) {
  # Network generation only
  # Get list of available projects with WGCNA results
  project_dirs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
  
  if (length(project_dirs) == 0) {
    cat("No projects with WGCNA results found. Please run WGCNA analysis first.\n")
  } else {
    cat("Found", length(project_dirs), "projects with WGCNA results:\n")
    print(project_dirs)
    
    selected_project_id <- readline(prompt = "Enter a project ID to generate network (e.g., TCGA-THYM): ")
    if (selected_project_id %in% project_dirs) {
      # Load WGCNA results
      module_file <- file.path(output_dir, selected_project_id, paste0(selected_project_id, "_gene_modules.csv"))
      if (!file.exists(module_file)) {
        cat("Module assignment file not found for", selected_project_id, "\n")
      } else {
        # Run only the network generation part
        # This is a bit hacky but we'll just call our full function with a flag
        result <- run_wgcna_workflow(selected_project_id, generate_network = TRUE)
      }
    } else {
      cat("Project", selected_project_id, "not found. Available projects:", paste(project_dirs, collapse=", "), "\n")
    }
  }
} else if (choice == 3) {
  # Full workflow
  selected_project_id <- readline(prompt = "Enter a project ID to analyze (e.g., TCGA-THYM): ")
  if (selected_project_id %in% project_ids) {
    result <- run_wgcna_workflow(selected_project_id, generate_network = TRUE)
    cat("\nComplete workflow finished for", selected_project_id, "\n")
  } else {
    cat("Project", selected_project_id, "not found. Available projects:", paste(project_ids, collapse=", "), "\n")
  }
} else {
  cat("Invalid choice. Please restart the script and choose 1, 2, or 3.\n")
}