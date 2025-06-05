# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidygraph)
library(stringr)
library(patchwork)

# Set working directory
setwd("/media/windows/BIOINFORMATICS/TCGA/TCGA-analysis/TCGA-Chaperones/heatmap_patients_stratification")

# Read the chaperone gene list for highlighting
gene_list <- read_csv("../gene_list.csv", show_col_types = FALSE)
chaperone_genes <- gene_list$Name

# Function to create tidygraph object from pathfindR data
create_pathway_graph <- function(pathway_data, limit_pathways = 8) {
  # Limit to top N pathways
  pathway_data <- head(pathway_data, limit_pathways)
  
  # Initialize edges and nodes dataframes
  edges <- data.frame(
    from = character(),
    to = character(),
    direction = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each pathway
  for (i in 1:nrow(pathway_data)) {
    pathway_term <- pathway_data$Term_Description[i]
    pathway_id <- pathway_data$ID[i]
    
    # Process up-regulated genes
    if (!is.na(pathway_data$Up_regulated[i])) {
      up_genes <- strsplit(pathway_data$Up_regulated[i], ",")[[1]] %>% trimws()
      
      if (length(up_genes) > 0) {
        up_edges <- data.frame(
          from = rep(pathway_term, length(up_genes)),
          to = up_genes,
          direction = "up",
          stringsAsFactors = FALSE
        )
        edges <- rbind(edges, up_edges)
      }
    }
    
    # Process down-regulated genes
    if (!is.na(pathway_data$Down_regulated[i])) {
      down_genes <- strsplit(pathway_data$Down_regulated[i], ",")[[1]] %>% trimws()
      
      if (length(down_genes) > 0) {
        down_edges <- data.frame(
          from = rep(pathway_term, length(down_genes)),
          to = down_genes,
          direction = "down",
          stringsAsFactors = FALSE
        )
        edges <- rbind(edges, down_edges)
      }
    }
  }
  
  # Create nodes data
  all_pathway_terms <- unique(edges$from)
  all_genes <- unique(edges$to)
  
  pathway_nodes <- data.frame(
    name = all_pathway_terms,
    type = "pathway",
    is_chaperone = FALSE,
    regulation = NA_character_,
    stringsAsFactors = FALSE
  )
  
  gene_nodes <- data.frame(
    name = all_genes,
    type = "gene",
    is_chaperone = all_genes %in% chaperone_genes,
    regulation = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Update gene regulation (some genes might be both up and down in different pathways)
  for (gene in all_genes) {
    current_index <- which(gene_nodes$name == gene)
    
    # Check if gene appears in up-regulated list
    if (gene %in% edges$to[edges$direction == "up"]) {
      gene_nodes$regulation[current_index] <- "up"
    }
    
    # Check if gene appears in down-regulated list
    if (gene %in% edges$to[edges$direction == "down"]) {
      # Fix the error: properly check for NA or "up" values
      current_regulation <- gene_nodes$regulation[current_index]
      
      # If already marked as up, it's in both categories
      if (!is.na(current_regulation) && current_regulation == "up") {
        gene_nodes$regulation[current_index] <- "both"
      } else {
        gene_nodes$regulation[current_index] <- "down"
      }
    }
  }
  
  # Combine all nodes
  all_nodes <- rbind(pathway_nodes, gene_nodes)
  
  # Create tidygraph object
  graph <- tbl_graph(
    nodes = all_nodes,
    edges = select(edges, from, to),
    directed = FALSE
  )
  
  # Add edge regulation type to graph edges
  edge_regulation <- edges$direction
  graph <- graph %>%
    activate(edges) %>%
    mutate(regulation = edge_regulation)
  
  return(graph)
}

# Function to create highly readable pathway graphs
create_readable_pathway_plot <- function(graph, project_name) {
  # Extract chaperone genes for later reference
  chaperone_genes_in_graph <- graph %>%
    activate(nodes) %>%
    filter(is_chaperone) %>%
    pull(name)
  
  # Create ggraph plot with careful layout
  pathway_plot <- graph %>%
    # Use nicely layout for better spacing
    ggraph(layout = "fr", niter = 10000) + 
    # Add edges with minimal alpha to reduce visual clutter
    geom_edge_link(
      aes(color = regulation),
      width = 0.6, 
      alpha = 0.7,
      show.legend = TRUE
    ) +
    # Add nodes with custom sizes and colors
    geom_node_point(
      aes(
        size = type,
        color = case_when(
          is_chaperone ~ "chaperone",
          type == "pathway" ~ "pathway",
          TRUE ~ regulation
        ),
        shape = type
      ),
      show.legend = TRUE
    ) +
    # Add labels with careful positioning - using fixed size directly instead of aes()
    geom_node_text(
      aes(
        label = name,
        fontface = ifelse(is_chaperone, "bold", "plain")
      ),
      size = 5,  # Use fixed size instead of mapping
      repel = TRUE,
      point.padding = unit(1.2, "lines"),
      box.padding = unit(1.5, "lines"),
      force = 10,
      max.overlaps = Inf,
      check_overlap = TRUE
    ) +
    # Define custom scales
    scale_shape_manual(
      values = c(pathway = 19, gene = 19),  # Changed gene shape to circle (19) instead of triangle (17)
      name = "Node Type"
    ) +
    scale_size_manual(
      values = c(pathway = 18, gene = 10),
      name = "Node Type"
    ) +
    scale_color_manual(
      values = c(
        "chaperone" = "#9932CC",   # Changed to purple for chaperones
        "pathway" = "#4682B4",     # Steel blue for pathways
        "up" = "#32CD32",          # Green for up-regulated genes
        "down" = "#E41A1C",        # Red for down-regulated genes (changed from purple)
        "both" = "#FFD700"         # Gold for genes regulated both ways
      ),
      name = "Node Type",
      labels = c(
        "chaperone" = "Chaperone Gene",
        "pathway" = "Pathway",
        "up" = "Up-regulated Gene",
        "down" = "Down-regulated Gene",
        "both" = "Bi-directional Gene"
      ),
      na.value = "gray50"
    ) +
    # Edge color scale
    scale_edge_color_manual(
      values = c(
        "up" = "#32CD32",    # Green for upregulated connections
        "down" = "#E41A1C"   # Red for downregulated connections
      ),
      name = "Regulation"
    ) +
    # Additional theme elements
    labs(
      title = paste("Chaperone-Related Pathways:", project_name),
      subtitle = paste0("Chaperone genes: ", paste(chaperone_genes_in_graph, collapse = ", ")),
      caption = "Purple: chaperone genes, Blue: pathways, Green: up-regulated, Red: down-regulated"
    ) +
    theme_graph(base_family = "sans", base_size = 16) +
    theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 20, hjust = 0.5),
      plot.margin = margin(50, 50, 50, 50),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold")
    )
  
  return(pathway_plot)
}

# Function to create a simplified plot focusing only on chaperones
create_chaperone_focused_plot <- function(graph, project_name) {
  # Extract chaperone nodes and connected pathway nodes
  chaperone_nodes <- graph %>%
    activate(nodes) %>%
    filter(is_chaperone) %>%
    pull(name)
  
  # If no chaperones, return NULL
  if (length(chaperone_nodes) == 0) {
    return(NULL)
  }
  
  # Get direct pathway connections to chaperones
  chap_connections <- graph %>%
    activate(edges) %>%
    filter(.N()$is_chaperone[to] | .N()$is_chaperone[from]) %>%
    as_tibble()
  
  # Create a subgraph with only chaperones and connected pathways
  subgraph_nodes <- unique(c(chaperone_nodes, 
                           chap_connections$from[!chap_connections$from %in% chaperone_nodes],
                           chap_connections$to[!chap_connections$to %in% chaperone_nodes]))
  
  focused_graph <- graph %>%
    activate(nodes) %>%
    filter(name %in% subgraph_nodes) %>%
    # Keep only relevant edges
    activate(edges) %>%
    filter(from %in% subgraph_nodes & to %in% subgraph_nodes)
  
  # Create simplified plot with a manual layout
  set.seed(42) # For reproducible layout
  
  # Create a text size variable in graph for use with filter
  focused_graph <- focused_graph %>%
    activate(nodes) %>%
    mutate(
      node_text_size = case_when(
        is_chaperone ~ 8,
        type == "pathway" ~ 10,
        TRUE ~ 0
      ),
      text_to_display = ifelse(node_text_size > 0, name, "")
    )
  
  chap_plot <- focused_graph %>%
    ggraph(layout = "kk") +
    geom_edge_link(
      aes(color = regulation),
      width = 1.2, 
      alpha = 0.7
    ) +
    geom_node_point(
      aes(
        size = type,
        color = case_when(
          is_chaperone ~ "chaperone",
          type == "pathway" ~ "pathway",
          TRUE ~ "other"
        ),
        shape = type
      ),
      show.legend = TRUE
    ) +
    # Only show labels for chaperones and pathways
    geom_node_text(
      aes(
        label = text_to_display,
        fontface = ifelse(is_chaperone, "bold", "plain")
      ),
      size = 6,  # Fixed size for all labels
      repel = TRUE, 
      point.padding = unit(2, "lines"),
      box.padding = unit(2, "lines"),
      force = 15
    ) +
    scale_shape_manual(
      values = c(pathway = 19, gene = 19),  # Changed to circles for genes as well
      name = "Type"
    ) +
    scale_size_manual(
      values = c(pathway = 20, gene = 15),
      name = "Type"
    ) +
    scale_color_manual(
      values = c(
        "chaperone" = "#9932CC",   # Changed to purple for chaperones
        "pathway" = "#4682B4",     # Blue for pathways
        "other" = "#A9A9A9"        # Gray for other genes
      ),
      name = "Type"
    ) +
    # Edge color scale
    scale_edge_color_manual(
      values = c(
        "up" = "#32CD32",    # Green for upregulated connections
        "down" = "#E41A1C"   # Red for downregulated connections
      ),
      name = "Regulation"
    ) +
    labs(
      title = paste("Chaperone-Pathway Interactions:", project_name),
      subtitle = paste("Focused on chaperone genes and their associated pathways")
    ) +
    theme_graph(base_family = "sans", base_size = 18) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 22, hjust = 0.5),
      plot.margin = margin(60, 60, 60, 60),
      plot.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20, face = "bold")
    )
  
  return(chap_plot)
}

# Main function to process all files
process_pathfindR_results <- function() {
  # Set paths
  filtered_dir <- "pathfindR_results/chaperone_filtered"
  output_dir <- "pathfindR_results/high_quality_graphs"
  
  # Check for filtered directory
  if (!dir.exists(filtered_dir)) {
    cat("Error: Filtered directory not found. Please run the filtering script first.\n")
    return(FALSE)
  }
  
  # List all filtered CSV files
  csv_files <- list.files(filtered_dir, pattern = "_chaperone_filtered\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cat("No filtered CSV files found in", filtered_dir, "\n")
    return(FALSE)
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Process each file
  for (csv_file in csv_files) {
    # Extract project name
    project_name <- basename(csv_file) %>% str_remove("_chaperone_filtered\\.csv$")
    
    cat("\n-----------------------------\n")
    cat("Processing:", project_name, "\n")
    
    # Create project directory
    project_dir <- file.path(output_dir, project_name)
    if (!dir.exists(project_dir)) {
      dir.create(project_dir, recursive = TRUE)
    }
    
    # Read data
    pathway_data <- read_csv(csv_file, show_col_types = FALSE)
    
    if (nrow(pathway_data) == 0) {
      cat("No data found for", project_name, "\n")
      next
    }
    
    # Create graph object
    graph <- create_pathway_graph(pathway_data, limit_pathways = 6) # Reduced from 8 to 6 for better readability
    
    # Generate the full network plot
    cat("Creating full network plot...\n")
    full_plot <- create_readable_pathway_plot(graph, project_name)
    
    # Save full plot
    ggsave(
      file.path(project_dir, paste0(project_name, "_full_network.pdf")),
      plot = full_plot,
      width = 30,
      height = 24,
      device = "pdf",
      limitsize = FALSE
    )
    
    ggsave(
      file.path(project_dir, paste0(project_name, "_full_network.png")),
      plot = full_plot,
      width = 30,
      height = 24,
      dpi = 300,
      device = "png",
      limitsize = FALSE
    )
    
    # Generate the chaperone-focused plot
    cat("Creating chaperone-focused plot...\n")
    chap_plot <- create_chaperone_focused_plot(graph, project_name)
    
    if (!is.null(chap_plot)) {
      # Save chaperone-focused plot
      ggsave(
        file.path(project_dir, paste0(project_name, "_chaperone_focused.pdf")),
        plot = chap_plot,
        width = 30,
        height = 24,
        device = "pdf",
        limitsize = FALSE
      )
      
      ggsave(
        file.path(project_dir, paste0(project_name, "_chaperone_focused.png")),
        plot = chap_plot,
        width = 30,
        height = 24,
        dpi = 300,
        device = "png",
        limitsize = FALSE
      )
    }
    
    cat("Completed processing", project_name, "\n")
  }
  
  cat("\nAll plots generated successfully!\n")
  cat("Output saved to:", output_dir, "\n")
  return(TRUE)
}

# Run the script
cat("Starting high-quality pathway visualization...\n")
cat("Creating graphs with MUCH better spacing and readability...\n\n")

process_pathfindR_results()

cat("\nVisualization complete!\n")
cat("Two types of plots were created for each project:\n")
cat("1. Full network with all genes and pathways\n")
cat("2. Chaperone-focused network showing only chaperones and their pathways\n")
cat("With updated colors: green (upregulated), red (downregulated), purple (chaperones)\n")