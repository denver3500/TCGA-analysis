# Load required packages
library(pathfindR)
library(dplyr)
library(readr)
library(ggkegg)
library(ggplot2)
library(stringr)

# Define input and output directories
deg_dir <- "TCGA-Chaperones/heatmap_patients_stratification/DEG_results"
output_dir <- "TCGA-Chaperones/heatmap_patients_stratification/pathfindR_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Find all DEG files - use new pattern for cluster comparisons
deg_files <- list.files(deg_dir, pattern = "cluster.*_vs_cluster.*_DEG\\.csv$", full.names = TRUE)

if (length(deg_files) == 0) {
  stop("No cluster comparison DEG files found. Make sure the files match the pattern 'cluster*_vs_cluster*_DEG.csv'")
}

# Extract project IDs and cluster comparison info
file_info <- str_match(deg_files, ".*/((TCGA-[A-Z]+)_cluster([0-9]+)_vs_cluster([0-9]+))_DEG\\.csv$")
comparison_ids <- file_info[, 2]  # Full comparison ID (e.g., "TCGA-BRCA_cluster2_vs_cluster1")
project_ids <- file_info[, 3]     # Just project ID (e.g., "TCGA-BRCA")
cluster1_ids <- file_info[, 5]    # Second cluster number
cluster2_ids <- file_info[, 4]    # First cluster number (reversed because file names are ordered cluster2_vs_cluster1)

cat("Found", length(deg_files), "DEG files to process\n")
cat("Projects represented:", length(unique(project_ids)), "\n")

# Create a summary dataframe
summary_results <- list()

# Process each DEG file
for (i in 1:length(deg_files)) {
  comparison_id <- comparison_ids[i]
  project_id <- project_ids[i]
  cluster1 <- cluster1_ids[i]
  cluster2 <- cluster2_ids[i]
  deg_file <- deg_files[i]
  
  cat("\n===== Processing", project_id, "| Cluster", cluster2, "vs Cluster", cluster1, "=====\n")
  start_time <- Sys.time()
  
  # Create project/comparison folder
  comparison_dir <- file.path(output_dir, comparison_id)
  dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Read and process DEG data
  deg_data <- read_csv(deg_file, show_col_types = FALSE)
  cat("Reading file:", basename(deg_file), "\n")
  
  # Filter for significant DEGs
  filtered_data <- deg_data %>% filter(FDR < 0.05)
  sig_count <- nrow(filtered_data)
  cat("Found", sig_count, "significant DEGs\n")
  
  # Skip if not enough DEGs
  if (sig_count < 10) {
    cat("Not enough significant DEGs (minimum 10). Skipping.\n")
    summary_results[[i]] <- list(
      Project = project_id,
      Comparison = paste0("Cluster", cluster2, " vs Cluster", cluster1),
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      Significant_DEGs = sig_count,
      Enriched_Terms = 0,
      Top_Term = "Insufficient DEGs",
      Top_Term_PValue = NA,
      Runtime_Minutes = 0
    )
    next
  }
  
  # Check if gene_name column exists
  if (!"gene_name" %in% colnames(filtered_data)) {
    cat("Error: gene_name column not found in DEG file. Available columns:", paste(colnames(filtered_data), collapse = ", "), "\n")
    summary_results[[i]] <- list(
      Project = project_id,
      Comparison = paste0("Cluster", cluster2, " vs Cluster", cluster1),
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      Significant_DEGs = sig_count,
      Enriched_Terms = 0,
      Top_Term = "Error: gene_name column not found",
      Top_Term_PValue = NA,
      Runtime_Minutes = 0
    )
    next
  }
  
  # Handle missing gene names
  missing_names <- sum(is.na(filtered_data$gene_name))
  if (missing_names > 0) {
    cat("Warning:", missing_names, "genes have missing gene names and will be excluded\n")
    filtered_data <- filtered_data %>% filter(!is.na(gene_name))
    sig_count <- nrow(filtered_data)
    
    # Check if we still have enough genes
    if (sig_count < 10) {
      cat("Not enough genes with valid names (minimum 10). Skipping.\n")
      summary_results[[i]] <- list(
        Project = project_id,
        Comparison = paste0("Cluster", cluster2, " vs Cluster", cluster1),
        Cluster1 = cluster1,
        Cluster2 = cluster2,
        Significant_DEGs = sig_count,
        Enriched_Terms = 0,
        Top_Term = "Insufficient genes with valid names",
        Top_Term_PValue = NA,
        Runtime_Minutes = 0
      )
      next
    }
  }
  
  # Create input for pathfindR
  input_df <- data.frame(
    Gene_symbol = filtered_data$gene_name,
    logFC = filtered_data$logFC,
    FDR_p = filtered_data$FDR
  )
  
  # Run pathfindR
  cat("Running pathfindR analysis...\n")
  res <- tryCatch({
    pathfindR_results <- run_pathfindR(
      input_df,
      gene_sets = "GO-All",
      pin_name = "KEGG",
      min_gset_size = 10,
      max_gset_size = 300,
      p_val_threshold = 0.05,
      output_dir = comparison_dir
    )
    
    # Save results
    result_file <- file.path(comparison_dir, paste0(comparison_id, "_pathfindR_results.csv"))
    write_csv(pathfindR_results, result_file)
    
    # Get term count and top term
    term_count <- nrow(pathfindR_results)
    if (term_count > 0) {
      top_term <- pathfindR_results$Term_Description[1]
      p_value <- pathfindR_results$lowest_p[1]
      cat("Found", term_count, "enriched terms\n")
      cat("Top term:", top_term, "(p =", p_value, ")\n")
      
      # Generate visualizations
      cat("Clustering enriched terms...\n")
      cluster_enriched_terms(pathfindR_results, output_dir = comparison_dir)
      
      # Create enrichment chart
      cat("Creating enrichment chart for top terms...\n")
      if (term_count > 1) {
        # Only show top 15 terms if there are more
        terms_to_chart <- min(15, term_count)
        enrichment_chart <- pathfindR_results[1:terms_to_chart,] %>%
          ggplot(aes(x = reorder(Term_Description, -Fold_Enrichment), 
                     y = Fold_Enrichment, 
                     fill = -log10(lowest_p))) +
          geom_bar(stat = "identity") +
          scale_fill_continuous(name = "-log10(p-value)") +
          labs(
            title = paste(project_id, "- Cluster", cluster2, "vs Cluster", cluster1),
            subtitle = "Top Enriched Terms",
            x = "",
            y = "Fold Enrichment"
          ) +
          coord_flip() +
          theme_bw() +
          theme(
            axis.text.y = element_text(size = 9),
            plot.title = element_text(size = 11, face = "bold"),
            legend.position = "right"
          )
        
        ggsave(file.path(comparison_dir, paste0(comparison_id, "_top_terms.pdf")), 
               enrichment_chart, width = 10, height = 7)
      }
      
      list(
        Project = project_id,
        Comparison = paste0("Cluster", cluster2, " vs Cluster", cluster1),
        Cluster1 = cluster1,
        Cluster2 = cluster2,
        Significant_DEGs = sig_count,
        Enriched_Terms = term_count,
        Top_Term = top_term,
        Top_Term_PValue = p_value,
        Runtime_Minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
      )
    } else {
      cat("No enriched terms found\n")
      list(
        Project = project_id,
        Comparison = paste0("Cluster", cluster2, " vs Cluster", cluster1),
        Cluster1 = cluster1,
        Cluster2 = cluster2,
        Significant_DEGs = sig_count, 
        Enriched_Terms = 0,
        Top_Term = "No enriched terms",
        Top_Term_PValue = NA,
        Runtime_Minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
      )
    }
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    list(
      Project = project_id,
      Comparison = paste0("Cluster", cluster2, " vs Cluster", cluster1),
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      Significant_DEGs = sig_count,
      Enriched_Terms = 0,
      Top_Term = paste("Error:", substr(e$message, 1, 50)),
      Top_Term_PValue = NA,
      Runtime_Minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
    )
  })
  
  summary_results[[i]] <- res
}

# Convert results to dataframe
summary_df <- do.call(rbind.data.frame, lapply(summary_results, as.data.frame))

# Save summary
summary_file <- file.path(output_dir, "pathfindR_analysis_summary.csv")
write_csv(summary_df, summary_file)

# Generate summary heatmaps by project
if (nrow(summary_df) > 0 && sum(summary_df$Enriched_Terms > 0) > 0) {
  # Get project list
  unique_projects <- unique(summary_df$Project)
  
  for (proj in unique_projects) {
    proj_data <- summary_df %>% filter(Project == proj)
    
    # Skip projects with no enriched terms
    if (sum(proj_data$Enriched_Terms) == 0) next
    
    # Create directory for project
    proj_dir <- file.path(output_dir, proj)
    dir.create(proj_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Get all enrichment results for this project
    enrichment_files <- list()
    for (i in 1:nrow(proj_data)) {
      comp_id <- paste0(proj, "_cluster", proj_data$Cluster2[i], "_vs_cluster", proj_data$Cluster1[i])
      result_file <- file.path(output_dir, comp_id, paste0(comp_id, "_pathfindR_results.csv"))
      
      if (file.exists(result_file)) {
        enrichment_files[[comp_id]] <- result_file
      }
    }
    
    # Generate comparison heatmap if multiple comparisons
    if (length(enrichment_files) > 1) {
      cat("\nGenerating comparison heatmap for", proj, "...\n")
      
      # Try to generate comparison visualization
      tryCatch({
        combined_data <- lapply(names(enrichment_files), function(comp_name) {
          df <- read_csv(enrichment_files[[comp_name]], show_col_types = FALSE)
          if (nrow(df) > 0) {
            df$Comparison <- comp_name
            return(df)
          } else {
            return(NULL)
          }
        })
        
        combined_data <- bind_rows(combined_data)
        
        # Keep only top terms
        top_terms <- combined_data %>%
          group_by(Term_Description) %>%
          summarise(min_p = min(lowest_p)) %>%
          arrange(min_p) %>%
          slice_head(n = 30) %>%
          pull(Term_Description)
        
        # Filter and reshape for heatmap
        heatmap_data <- combined_data %>%
          filter(Term_Description %in% top_terms) %>%
          select(Term_Description, Comparison, lowest_p) %>%
          mutate(neg_log_p = -log10(lowest_p)) %>%
          tidyr::pivot_wider(
            id_cols = Term_Description,
            names_from = Comparison,
            values_from = neg_log_p,
            values_fill = 0
          ) %>%
          as.data.frame()
        
        # Set rownames for pheatmap
        rownames(heatmap_data) <- heatmap_data$Term_Description
        heatmap_data <- heatmap_data[, -1]
        
        # Generate heatmap
        pdf(file.path(proj_dir, paste0(proj, "_pathway_comparison.pdf")), width = 12, height = 14)
        pheatmap::pheatmap(
          heatmap_data,
          main = paste(proj, "- Pathway Enrichment Across Cluster Comparisons"),
          color = colorRampPalette(c("white", "red"))(100),
          fontsize_row = 8,
          fontsize_col = 8,
          angle_col = 45,
          clustering_method = "ward.D2"
        )
        dev.off()
      }, error = function(e) {
        cat("Error generating comparison heatmap for", proj, ":", e$message, "\n")
      })
    }
  }
}

# Print final summary
cat("\n===== PathfindR Analysis Summary =====\n")
cat("Projects processed:", length(unique(summary_df$Project)), "\n")
cat("Cluster comparisons processed:", nrow(summary_df), "\n")
cat("Comparisons with enriched terms:", sum(summary_df$Enriched_Terms > 0), "\n")
cat("Total enriched terms:", sum(summary_df$Enriched_Terms), "\n")
cat("Average runtime per comparison:", round(mean(summary_df$Runtime_Minutes, na.rm = TRUE), 2), "minutes\n")

# Print the summary table
print(summary_df)

cat("\nPathfindR analysis complete for all cluster comparisons!\n")