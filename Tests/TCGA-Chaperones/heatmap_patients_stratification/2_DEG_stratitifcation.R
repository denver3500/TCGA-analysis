library(readr)
library(dplyr)
library(SummarizedExperiment)
library(edgeR)
library(stringr)
library(tidyr)
library(ggplot2)
library(purrr)
library(RColorBrewer)
library(pheatmap)

# Define file paths
rds_dir <- "TCGA-Chaperones/rds"
cluster_dir <- "TCGA-Chaperones/heatmap_patients_stratification/deseq2_clustering" 
output_dir <- "TCGA-Chaperones/heatmap_patients_stratification/DEG_results"
log_file <- "TCGA-Chaperones/heatmap_patients_stratification/DEG_analysis_log.txt"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
con <- file(log_file, "w")
sink(con, split = TRUE)

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
}

# Function to perform DEG analysis for a project using clusters
analyze_project_clusters <- function(project_id) {
  message("\n====== Processing ", project_id, " ======")
  
  # File paths for cluster assignments
  cluster_file <- file.path(cluster_dir, project_id, paste0(project_id, "_patient_clusters.csv"))
  
  # Check if cluster file exists
  if (!file.exists(cluster_file)) {
    message("  Cluster file not found for ", project_id, ". Skipping.")
    return(NULL)
  }
  
  # Read cluster assignments
  message("  Loading cluster assignments")
  cluster_data <- read_csv(cluster_file, show_col_types = FALSE)
  
  # Check and ensure Patient_Cluster is properly formatted
  if ("Patient_Cluster" %in% colnames(cluster_data)) {
    # Convert to character if needed
    if (is.list(cluster_data$Patient_Cluster)) {
      message("  Converting Patient_Cluster from list to character vector")
      cluster_data$Patient_Cluster <- as.character(unlist(cluster_data$Patient_Cluster))
    } else {
      cluster_data$Patient_Cluster <- as.character(cluster_data$Patient_Cluster)
    }
  } else {
    message("  Error: Patient_Cluster column not found in cluster data")
    return(NULL)
  }
  
  # Get unique clusters
  unique_clusters <- sort(unique(cluster_data$Patient_Cluster))
  n_clusters <- length(unique_clusters)
  
  message("  Found ", n_clusters, " patient clusters")
  
  # Show cluster distribution using base R for robustness
  cluster_counts <- table(cluster_data$Patient_Cluster)
  message("  Cluster distribution:")
  print(cluster_counts)
  
  # Find and load the RDS file for this project
  rds_files <- list.files(rds_dir, pattern = paste0("^", project_id, "_.*\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("  No RDS file found for project ", project_id, ". Skipping.")
    return(NULL)
  }
  
  message("  Loading RDS file: ", basename(rds_files[1]))
  se_object <- readRDS(rds_files[1])
  
  # Check available assays in the SummarizedExperiment object
  message("  Available assays: ", paste(assayNames(se_object), collapse = ", "))
  
  # Choose the appropriate assay for expression data - prioritizing counts for DEG
  if ("counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "counts")
    message("  Using counts for differential expression analysis")
  } else if ("HTSeq - Counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "HTSeq - Counts")
    message("  Using HTSeq - Counts for differential expression analysis")
  } else if ("raw_counts" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "raw_counts")
    message("  Using raw_counts for differential expression analysis")
  } else if ("unstranded" %in% assayNames(se_object)) {
    count_data <- assay(se_object, "unstranded")
    message("  Using unstranded counts for differential expression analysis")
  } else {
    message("  No count data found. Cannot perform differential expression analysis.")
    return(NULL)
  }
  
  # Get patient information
  tcga_barcodes <- colnames(count_data)
  short_barcodes <- substr(tcga_barcodes, 1, 15)  # Use the first 15 characters for matching
  
  # Create a mapping between patient IDs in cluster data and RDS data
  cluster_data$short_barcode <- substr(cluster_data$Patient_ID, 1, 15)
  
  # Generate all pairwise combinations of clusters
  cluster_pairs <- combn(unique_clusters, 2, simplify = FALSE)
  message("  Will perform ", length(cluster_pairs), " pairwise cluster comparisons")
  
  # Initialize results list
  pairwise_results <- list()
  
  # Perform pairwise DEG analysis
  for (pair_idx in seq_along(cluster_pairs)) {
    cluster1 <- cluster_pairs[[pair_idx]][1]
    cluster2 <- cluster_pairs[[pair_idx]][2]
    
    message("\n  === Analyzing Cluster ", cluster1, " vs Cluster ", cluster2, " ===")
    
    # Get patients in each cluster
    cluster1_patients <- cluster_data$short_barcode[cluster_data$Patient_Cluster == cluster1]
    cluster2_patients <- cluster_data$short_barcode[cluster_data$Patient_Cluster == cluster2]
    
    # Match patients to RDS data indices
    cluster1_indices <- which(short_barcodes %in% cluster1_patients)
    cluster2_indices <- which(short_barcodes %in% cluster2_patients)
    
    # Check if we have enough samples
    if (length(cluster1_indices) < 3 || length(cluster2_indices) < 3) {
      message("  Not enough samples matched for clusters ", cluster1, " and ", cluster2, ". Skipping.")
      message("  Matched samples: ", length(cluster1_indices), " in cluster ", cluster1, 
              ", ", length(cluster2_indices), " in cluster ", cluster2)
      next
    }
    
    message("  Matched ", length(cluster1_indices), " samples in cluster ", cluster1, 
            " and ", length(cluster2_indices), " samples in cluster ", cluster2)
    
    # Create DGEList object
    sample_indices <- c(cluster1_indices, cluster2_indices)
    group <- factor(c(rep(paste0("Cluster", cluster1), length(cluster1_indices)), 
                     rep(paste0("Cluster", cluster2), length(cluster2_indices))), 
                   levels = c(paste0("Cluster", cluster1), paste0("Cluster", cluster2)))
    
    # Create a DGEList object
    dge <- DGEList(counts = count_data[, sample_indices], group = group)
    
    # Filter low expressed genes
    message("  Filtering lowly expressed genes")
    keep <- filterByExpr(dge)
    message("  Keeping ", sum(keep), " out of ", nrow(dge), " genes")
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    # Normalize the data
    message("  Normalizing data")
    dge <- calcNormFactors(dge)
    
    # Design matrix
    design <- model.matrix(~group)
    colnames(design) <- c("Intercept", paste0("Cluster", cluster2, "vsCluster", cluster1))
    
    # Estimate dispersion
    message("  Estimating dispersion")
    dge <- estimateDisp(dge, design)
    
    # Fit the model
    message("  Fitting GLM")
    fit <- glmQLFit(dge, design)
    
    # Test for differential expression
    message("  Testing for differential expression")
    qlf <- glmQLFTest(fit, coef = 2)  # Second cluster vs first cluster
    
    # Extract results
    message("  Extracting DEG results")
    res <- topTags(qlf, n = Inf)
    results_df <- as.data.frame(res)
    results_df$gene_ensembl <- rownames(results_df)
    
    # Add gene names if available in the SummarizedExperiment object
    if (is(se_object, "SummarizedExperiment") && is(try(rowRanges(se_object), silent = TRUE), "GRanges")) {
      gene_info <- as.data.frame(mcols(rowRanges(se_object)))
      
      # Find column with gene names
      potential_name_cols <- c("gene_name", "external_gene_name", "symbol", "gene_symbol")
      name_col <- potential_name_cols[potential_name_cols %in% colnames(gene_info)][1]
      
      if (!is.na(name_col)) {
        # Create mapping from ENSEMBL ID to gene name
        ensembl_to_gene_map <- data.frame(
          ensembl_id = rownames(se_object),
          gene_name = gene_info[[name_col]],
          stringsAsFactors = FALSE
        )
        
        # Extract base ENSEMBL IDs
        ensembl_to_gene_map$ensembl_base <- extract_base_ensembl(ensembl_to_gene_map$ensembl_id)
        results_df$ensembl_base <- extract_base_ensembl(results_df$gene_ensembl)
        
        # Merge gene names to results
        results_df <- merge(results_df, 
                           ensembl_to_gene_map[, c("ensembl_base", "gene_name")], 
                           by = "ensembl_base", all.x = TRUE)
      }
    }
    
    # If gene names not available, add an ensembl_base column
    if (!"gene_name" %in% colnames(results_df) && !"ensembl_base" %in% colnames(results_df)) {
      results_df$ensembl_base <- extract_base_ensembl(results_df$gene_ensembl)
    }
    
    # Reorder columns to put gene names first if they exist
    if ("gene_name" %in% colnames(results_df)) {
      results_df <- results_df %>%
        select(ensembl_base, gene_name, gene_ensembl, everything()) %>%
        arrange(PValue)
    } else {
      results_df <- results_df %>%
        select(ensembl_base, gene_ensembl, everything()) %>%
        arrange(PValue)
    }
    
    # Add comparison information
    results_df$comparison <- paste0("Cluster", cluster2, "_vs_Cluster", cluster1)
    results_df$cluster1 <- cluster1
    results_df$cluster2 <- cluster2
    
    # Save results
    comparison_name <- paste0("cluster", cluster2, "_vs_cluster", cluster1)
    results_file <- file.path(output_dir, paste0(project_id, "_", comparison_name, "_DEG.csv"))
    write_csv(results_df, results_file)
    
    # Create summary of results
    up_genes <- sum(results_df$logFC > 0 & results_df$FDR < 0.05)
    down_genes <- sum(results_df$logFC < 0 & results_df$FDR < 0.05)
    total_sig <- sum(results_df$FDR < 0.05)
    
    message("  DEG analysis complete. Results saved to ", results_file)
    message("  Summary: ", total_sig, " significant genes (FDR < 0.05)")
    message("           ", up_genes, " up-regulated in cluster ", cluster2)
    message("           ", down_genes, " down-regulated in cluster ", cluster2)
    
    # Create MA plot if there are significant DEGs
    if (total_sig > 0) {
      plot_file <- file.path(output_dir, paste0(project_id, "_", comparison_name, "_MAplot.pdf"))
      pdf(plot_file, width = 10, height = 8)
      
      # Create data frame for plotting
      plot_data <- results_df %>%
        mutate(significance = case_when(
          FDR < 0.01 ~ "FDR < 0.01",
          FDR < 0.05 ~ "FDR < 0.05",
          TRUE ~ "Not significant"
        ),
        significance = factor(significance, levels = c("FDR < 0.01", "FDR < 0.05", "Not significant")))
      
      # Create MA plot
      p <- ggplot(plot_data, aes(x = logCPM, y = logFC, color = significance)) +
        geom_point(size = 0.8, alpha = 0.7) +
        scale_color_manual(values = c("FDR < 0.01" = "red", "FDR < 0.05" = "orange", "Not significant" = "grey")) +
        labs(title = paste(project_id, "- Cluster", cluster2, "vs Cluster", cluster1),
             subtitle = paste("Up:", up_genes, "Down:", down_genes, "Total Sig:", total_sig),
             x = "Average Expression (logCPM)",
             y = "Log Fold Change (Cluster 2 vs Cluster 1)") +
        theme_bw() +
        theme(legend.position = "bottom")
      
      print(p)
      dev.off()
      message("  MA plot saved to ", plot_file)
      
      # Also create a volcano plot
      plot_file <- file.path(output_dir, paste0(project_id, "_", comparison_name, "_volcano.pdf"))
      pdf(plot_file, width = 10, height = 8)
      
      # Create volcano plot
      p <- ggplot(plot_data, aes(x = logFC, y = -log10(FDR), color = significance)) +
        geom_point(size = 0.8, alpha = 0.7) +
        scale_color_manual(values = c("FDR < 0.01" = "red", "FDR < 0.05" = "orange", "Not significant" = "grey")) +
        labs(title = paste(project_id, "- Cluster", cluster2, "vs Cluster", cluster1),
             subtitle = paste("Up:", up_genes, "Down:", down_genes, "Total Sig:", total_sig),
             x = "Log Fold Change (Cluster 2 vs Cluster 1)",
             y = "-log10 FDR") +
        theme_bw() +
        theme(legend.position = "bottom") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey")
      
      print(p)
      dev.off()
      message("  Volcano plot saved to ", plot_file)
      
      # Store summary of this comparison
      pairwise_results[[paste0("Cluster", cluster2, "_vs_Cluster", cluster1)]] <- data.frame(
        Project = project_id,
        Comparison = paste0("Cluster", cluster2, "_vs_Cluster", cluster1),
        Cluster1 = cluster1,
        Cluster2 = cluster2,
        Cluster1_Samples = length(cluster1_indices),
        Cluster2_Samples = length(cluster2_indices),
        Filtered_Genes = sum(keep),
        Total_DEGs = total_sig,
        Up_in_Cluster2 = up_genes,
        Down_in_Cluster2 = down_genes,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all pairwise results
  if (length(pairwise_results) > 0) {
    combined_results <- bind_rows(pairwise_results)
    
    # Save project-level summary
    project_summary_file <- file.path(output_dir, paste0(project_id, "_DEG_summary.csv"))
    write_csv(combined_results, project_summary_file)
    message("\n  Project summary saved to ", project_summary_file)
    
    # Return combined results
    return(combined_results)
  } else {
    message("\n  No significant DEGs found in any cluster comparison")
    return(NULL)
  }
}

# Find all project directories in the clustering output
project_dirs <- list.dirs(cluster_dir, full.names = FALSE, recursive = FALSE)
project_dirs <- project_dirs[project_dirs != ""] # Remove empty entries

message("Found ", length(project_dirs), " projects with cluster data")
print(project_dirs)

# Process each project
results <- list()

for (project_id in project_dirs) {
  result <- try(analyze_project_clusters(project_id))
  if (!inherits(result, "try-error") && !is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Generate summary table
if (length(results) > 0) {
  summary_df <- bind_rows(results)
  write_csv(summary_df, file.path(output_dir, "DEG_analysis_summary.csv"))
  
  message("\n====== Overall Summary ======")
  message("Total projects processed: ", length(unique(summary_df$Project)))
  message("Total pairwise comparisons: ", nrow(summary_df))
  message("Total comparisons with DEGs: ", sum(summary_df$Total_DEGs > 0))
  message("Total significant DEGs across all comparisons: ", sum(summary_df$Total_DEGs))
  message("Average DEGs per comparison: ", round(mean(summary_df$Total_DEGs), 1))
  
  # Create summary heatmap of DEG counts
  heatmap_data <- summary_df %>%
    select(Project, Comparison, Total_DEGs) %>%
    tidyr::pivot_wider(names_from = Comparison, values_from = Total_DEGs, values_fill = 0) %>%
    as.data.frame()
  
  rownames(heatmap_data) <- heatmap_data$Project
  heatmap_data <- heatmap_data[, -1]
  
  if (ncol(heatmap_data) > 0 && nrow(heatmap_data) > 0) {
    pdf(file.path(output_dir, "DEG_counts_heatmap.pdf"), width = 12, height = 10)
    pheatmap::pheatmap(
      heatmap_data,
      main = "Number of DEGs by Cluster Comparison",
      color = colorRampPalette(c("white", "red"))(100),
      display_numbers = TRUE,
      fontsize_number = 8,
      fontsize = 10,
      fontsize_row = 10,
      fontsize_col = 8,
      angle_col = 45
    )
    dev.off()
    message("Generated heatmap of DEG counts")
  }
  
  # Print summary table
  print(summary_df)
} else {
  message("No projects were successfully processed")
}

# Close log connection
sink()
close(con)

message("\nDifferential expression analysis based on patient clusters complete!")