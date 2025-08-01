# Dendrogram Analysis for Pan-Cancer Patient Stratification
# This script focuses on generating dendrograms and heatmaps with tree cuts
# Uses existing caching logic and saves images to pictures directory

# Load required libraries
library(tibble)
library(readr)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dendextend)  # For dendrogram manipulation
library(cluster)     # For clustering evaluation
library(viridis)     # For colorblind-friendly palettes
library(stringr)
library(circlize)    # For colorRamp2 function

# Define file paths - using paths from project root
gene_list_file <- "TCGA-Chaperones/gene_list.csv"
raw_data_dir <- "raw_data"
output_dir <- "TCGA-Chaperones/Patient_stratification/pan_cancer_clustering"
pictures_dir <- "TCGA-Chaperones/Patient_stratification/pictures"

# Create directories if they don't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pictures_dir, showWarnings = FALSE, recursive = TRUE)

# Define cache files (using existing logic)
cache_dir <- file.path(output_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

combined_data_cache <- file.path(cache_dir, "combined_expression_data.rds")
clustering_cache <- file.path(cache_dir, "clustering_results.rds")

message("=== Dendrogram Analysis for Pan-Cancer Patient Stratification ===")
message("Starting at: ", Sys.time())

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
}

# Function to process a single project and extract chaperone expression
process_single_project <- function(project_id, gene_list) {
  message("  Processing project: ", project_id)
  
  # Find and load RDS file from raw_data directory
  rds_files <- list.files(raw_data_dir, pattern = paste0("^", project_id, "_.*transcriptomic_exp\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("    No RDS file found for project ", project_id)
    return(NULL)
  }
  
  se_object <- readRDS(rds_files[1])
  
  # Filter for tumor samples only
  if ("definition" %in% colnames(colData(se_object))) {
    se_object <- se_object[, se_object$definition == "Primary solid Tumor"]
  } else if ("sample_type" %in% colnames(colData(se_object))) {
    se_object <- se_object[, se_object$sample_type == "Primary Tumor"]
  } else if (any(grepl("shortLetterCode", colnames(colData(se_object))))) {
    se_object <- se_object[, se_object$shortLetterCode == "TP"]
  }
  
  n_samples <- ncol(se_object)
  message("    Tumor samples: ", n_samples)
  
  if (n_samples < 10) {
    message("    Too few samples - skipping")
    return(NULL)
  }
  
  # Get RAW COUNT data for DESeq2
  assay_names <- assayNames(se_object)
  
  if ("unstranded" %in% assay_names) {
    count_data <- assay(se_object, "unstranded")
    used_assay <- "unstranded"
  } else if ("stranded_first" %in% assay_names) {
    count_data <- assay(se_object, "stranded_first")
    used_assay <- "stranded_first"
  } else if ("stranded_second" %in% assay_names) {
    count_data <- assay(se_object, "stranded_second")
    used_assay <- "stranded_second"
  } else {
    count_data <- assay(se_object, assay_names[1])
    used_assay <- assay_names[1]
    message("    WARNING: Using ", assay_names[1], " (may not be raw counts)")
  }
  
  # Round counts to integers
  count_data <- round(count_data)
  
  # Create sample metadata
  sample_metadata <- data.frame(
    sample_id = colnames(count_data),
    project = project_id,
    sample_type = "tumor",
    stringsAsFactors = FALSE
  )
  rownames(sample_metadata) <- colnames(count_data)
  
  # Create DESeq2 dataset with ALL GENES
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_metadata,
    design = ~ 1
  )
  
  # Filter low count genes
  min_samples <- ceiling(0.25 * ncol(dds))
  keep <- rowSums(counts(dds) >= 10) >= min_samples
  dds <- dds[keep, ]
  
  message("    Genes after filtering: ", nrow(dds))
  
  # DESeq2 normalization and VST transformation
  dds <- DESeq(dds, quiet = TRUE)
  vsd <- vst(dds, blind = TRUE)
  vst_all_genes <- assay(vsd)
  
  # Extract chaperone genes
  all_genes <- rownames(vst_all_genes)
  
  if (any(grepl("\\.", all_genes))) {
    all_base_genes <- extract_base_ensembl(all_genes)
    gene_mapping <- data.frame(
      full_id = all_genes,
      base_id = all_base_genes,
      stringsAsFactors = FALSE
    )
  } else {
    gene_mapping <- data.frame(
      full_id = all_genes,
      base_id = all_genes,
      stringsAsFactors = FALSE
    )
  }
  
  # Find chaperone genes
  chaperone_genes <- c()
  chaperone_names <- c()
  chaperone_families <- c()
  
  for (i in 1:nrow(gene_list)) {
    ensembl_id <- gene_list$ENSEMBL[i]
    gene_name <- gene_list$Name[i]
    family <- gene_list$Family[i]
    
    match_idx <- which(gene_mapping$base_id == ensembl_id)
    
    if (length(match_idx) > 0) {
      full_id <- gene_mapping$full_id[match_idx[1]]
      chaperone_genes <- c(chaperone_genes, full_id)
      chaperone_names <- c(chaperone_names, gene_name)
      chaperone_families <- c(chaperone_families, family)
    }
  }
  
  if (length(chaperone_genes) == 0) {
    message("    No chaperone genes found")
    return(NULL)
  }
  
  message("    Chaperone genes found: ", length(chaperone_genes))
  
  # Extract chaperone expression
  vst_chaperones <- vst_all_genes[chaperone_genes, , drop = FALSE]
  rownames(vst_chaperones) <- chaperone_names
  
  # Add project prefix to sample names
  colnames(vst_chaperones) <- paste0(project_id, ":", colnames(vst_chaperones))
  
  return(list(
    expression = vst_chaperones,
    project = project_id,
    n_samples = ncol(vst_chaperones),
    families = chaperone_families,
    used_assay = used_assay
  ))
}

# Function to combine all projects (using existing caching logic)
combine_all_projects <- function(force_recompute = FALSE) {
  
  # Check if cached data exists and is recent
  if (!force_recompute && file.exists(combined_data_cache)) {
    message("\n=== Loading cached combined expression data ===")
    combined_data <- readRDS(combined_data_cache)
    message("Loaded combined data with ", ncol(combined_data$combined_expression), " patients from ", length(combined_data$project_info), " projects")
    return(combined_data)
  }
  
  message("\n=== Loading and combining all TCGA projects ===")
  
  # Load gene list and get all TCGA .rds files from raw_data directory
  gene_list <- read_csv(gene_list_file, show_col_types = FALSE)
  tcga_files <- list.files(raw_data_dir, pattern = "^TCGA-.*_transcriptomic_exp\\.rds$", full.names = TRUE)
  
  message("Loaded gene list with ", nrow(gene_list), " genes")
  message("Found ", length(tcga_files), " TCGA files")
  
  # Extract project IDs from filenames
  extract_project_id <- function(filename) {
    basename(filename) %>%
      str_replace("_transcriptomic_exp\\.rds$", "")
  }
  
  project_ids <- sapply(tcga_files, extract_project_id)
  
  # Process all projects
  all_expression_data <- list()
  project_summaries <- list()
  
  for (i in seq_along(project_ids)) {
    project_id <- project_ids[i]
    
    tryCatch({
      project_result <- process_single_project(project_id, gene_list)
      if (!is.null(project_result)) {
        all_expression_data[[project_id]] <- project_result$expression
        project_summaries[[project_id]] <- list(
          project = project_id,
          n_samples = project_result$n_samples,
          used_assay = project_result$used_assay
        )
      }
    }, error = function(e) {
      message("    Error processing ", project_id, ": ", e$message)
    })
  }
  
  if (length(all_expression_data) == 0) {
    stop("No projects were successfully processed!")
  }
  
  message("\nSuccessfully processed ", length(all_expression_data), " projects")
  
  # Combine all expression matrices
  message("\n=== Combining expression matrices ===")
  
  # Get common genes across all projects
  all_gene_names <- Reduce(intersect, lapply(all_expression_data, rownames))
  message("Common chaperone genes across all projects: ", length(all_gene_names))
  
  if (length(all_gene_names) == 0) {
    stop("No common chaperone genes found across projects!")
  }
  
  # Combine matrices
  combined_matrices <- list()
  for (project_id in names(all_expression_data)) {
    expr_matrix <- all_expression_data[[project_id]]
    combined_matrices[[project_id]] <- expr_matrix[all_gene_names, , drop = FALSE]
  }
  
  # Merge all matrices
  combined_expression <- do.call(cbind, combined_matrices)
  
  message("Combined expression matrix: ", nrow(combined_expression), " genes × ", ncol(combined_expression), " patients")
  
  # Create patient metadata
  patient_metadata <- data.frame(
    Patient_ID = colnames(combined_expression),
    Project = str_extract(colnames(combined_expression), "^[^:]+"),
    Sample_ID = str_extract(colnames(combined_expression), "[^:]+$"),
    stringsAsFactors = FALSE
  )
  
  # Get gene families (use first project as reference)
  first_project_families <- gene_list$Family[match(all_gene_names, gene_list$Name)]
  names(first_project_families) <- all_gene_names
  
  # Combine everything
  combined_data <- list(
    combined_expression = combined_expression,
    patient_metadata = patient_metadata,
    gene_families = first_project_families,
    project_info = project_summaries,
    common_genes = all_gene_names
  )
  
  # Cache the combined data
  message("Caching combined data to: ", combined_data_cache)
  saveRDS(combined_data, combined_data_cache)
  
  return(combined_data)
}

# Function to perform clustering and generate dendrogram
perform_clustering <- function(combined_data, force_recompute = FALSE) {
  
  # Check if cached clustering exists and is compatible
  if (!force_recompute && file.exists(clustering_cache)) {
    message("\n=== Loading cached clustering results ===")
    tryCatch({
      clustering_results <- readRDS(clustering_cache)
      
      # Check if the cached results have the required clustering fields
      if (!is.null(clustering_results$patient_hclust) &&
          !is.null(clustering_results$patient_dend)) {
        
        # If expression matrix is missing, use it from combined_data
        if (is.null(clustering_results$expression_matrix)) {
          message("Expression matrix missing from cache, using combined_data expression matrix")
          clustering_results$expression_matrix <- combined_data$combined_expression
        }
        
        message("Loaded clustering results with expression matrix: ", 
                nrow(clustering_results$expression_matrix), " x ", 
                ncol(clustering_results$expression_matrix))
        return(clustering_results)
      } else {
        message("Cached clustering results are incomplete, recomputing...")
      }
    }, error = function(e) {
      message("Error loading cached clustering results: ", e$message)
      message("Recomputing clustering...")
    })
  }
  
  message("\n=== Performing patient clustering for dendrogram analysis ===")
  
  expression_matrix <- combined_data$combined_expression
  
  # Validate expression matrix
  if (is.null(expression_matrix) || nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    stop("Combined expression matrix is NULL or empty!")
  }
  
  message("Expression matrix for clustering: ", nrow(expression_matrix), " x ", ncol(expression_matrix))
  
  # Compute distance matrix for patients
  message("Computing patient distance matrix...")
  patient_dist <- dist(t(expression_matrix))
  
  # Perform hierarchical clustering
  message("Performing hierarchical clustering...")
  patient_hclust <- hclust(patient_dist, method = "ward.D2")
  
  # Create dendrogram
  patient_dend <- as.dendrogram(patient_hclust)
  
  # Store clustering results with expression matrix
  clustering_results <- list(
    patient_hclust = patient_hclust,
    patient_dend = patient_dend,
    patient_dist = patient_dist,
    expression_matrix = expression_matrix  # Ensure this is stored
  )
  
  # Cache clustering results
  message("Caching clustering results to: ", clustering_cache)
  saveRDS(clustering_results, clustering_cache)
  
  return(clustering_results)
}

# Function to generate dendrogram and heatmaps with 3 tree cuts
generate_dendrogram_analysis <- function(combined_data, clustering_results) {
  message("\n=== Generating dendrogram and heatmaps with tree cuts ===")
  
  expression_matrix <- clustering_results$expression_matrix
  patient_metadata <- combined_data$patient_metadata
  
  # Debug information
  message("Expression matrix dimensions: ", nrow(expression_matrix), " x ", ncol(expression_matrix))
  message("Patient metadata rows: ", nrow(patient_metadata))
  
  # Check for valid expression data
  if (is.null(expression_matrix) || nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    stop("Expression matrix is NULL or empty!")
  }
  
  # Check for NA or infinite values
  na_count <- sum(is.na(expression_matrix))
  inf_count <- sum(is.infinite(expression_matrix))
  
  if (na_count > 0 || inf_count > 0) {
    message("Warning: Expression matrix contains ", na_count, " NA values and ", inf_count, " infinite values")
    # Replace NA and infinite values with median
    median_val <- median(expression_matrix, na.rm = TRUE)
    expression_matrix[is.na(expression_matrix) | is.infinite(expression_matrix)] <- median_val
    message("Replaced NA/infinite values with median: ", median_val)
  }
  
  # Ensure expression matrix is numeric
  if (!is.numeric(expression_matrix)) {
    message("Converting expression matrix to numeric...")
    expression_matrix <- as.matrix(expression_matrix)
    mode(expression_matrix) <- "numeric"
  }
  
  # Define the 3 tree cuts
  tree_cuts <- c(3, 5, 7)
  
  # 1. Generate patient dendrogram
  message("Generating patient dendrogram...")
  
  dend_file <- file.path(pictures_dir, "patient_dendrogram.pdf")
  pdf(dend_file, width = 16, height = 8)
  
  # Plot dendrogram
  par(mar = c(6, 4, 4, 2))
  plot(clustering_results$patient_dend, 
       main = paste0("Pan-Cancer Patient Clustering Dendrogram (", 
                    ncol(expression_matrix), " patients)"), 
       horiz = FALSE, 
       axes = FALSE,
       cex = 0.6)
  
  dev.off()
  message("Generated dendrogram: ", dend_file)
  
  # 2. Generate heatmaps with each tree cut
  for (k in tree_cuts) {
    message("Generating heatmap with ", k, " clusters...")
    
    # Debug: Check expression matrix before each heatmap
    message("  Expression matrix for k=", k, ": ", nrow(expression_matrix), " x ", ncol(expression_matrix))
    message("  Class: ", class(expression_matrix), ", Mode: ", mode(expression_matrix))
    
    if (is.null(expression_matrix) || nrow(expression_matrix) == 0) {
      message("  Skipping k=", k, " - expression matrix is NULL or empty")
      next
    }
    
    # Cut tree to get clusters
    patient_clusters <- cutree(clustering_results$patient_hclust, k = k)
    
    # Create patient annotation (only clusters, no projects)
    patient_annotation <- data.frame(
      Cluster = factor(patient_clusters, levels = 1:k)
    )
    rownames(patient_annotation) <- patient_metadata$Patient_ID
    
    # Create colors (only for clusters)
    cluster_colors <- setNames(viridis(k), levels(patient_annotation$Cluster))
    
    annotation_colors <- list(
      Cluster = cluster_colors
    )
    
    # Create expression color scale using colorRamp2 for better control
    expression_range <- range(expression_matrix, na.rm = TRUE)
    expression_color_fun <- colorRamp2(
      c(expression_range[1], 
        expression_range[1] + (expression_range[2] - expression_range[1]) * 1/3,
        expression_range[1] + (expression_range[2] - expression_range[1]) * 2/3,
        expression_range[2]), 
      c("#fff0f3", "#ff8fa3", "#c9184a", "#590d22")
    )
    color_palette <- expression_color_fun(seq(expression_range[1], expression_range[2], length.out = 100))
    
    # Calculate gaps between clusters
    patient_order <- clustering_results$patient_hclust$order
    patient_cluster_ordered <- patient_clusters[patient_order]
    patient_cluster_gaps <- which(diff(patient_cluster_ordered) != 0)
    
    # Generate heatmap
    heatmap_file <- file.path(pictures_dir, paste0("heatmap_", k, "_clusters.pdf"))
    pdf(heatmap_file, width = 16, height = 8)  # Smaller PDF size
    
    # Show column names only for very small datasets
    show_colnames <- ncol(expression_matrix) <= 50
    
    tryCatch({
      # Ensure the matrix is properly formatted
      if (!is.matrix(expression_matrix)) {
        expression_matrix <- as.matrix(expression_matrix)
      }
      
      pheatmap(
        expression_matrix,
        color = color_palette,
        cluster_rows = TRUE,
        cluster_cols = clustering_results$patient_hclust,
        annotation_col = patient_annotation,
        annotation_colors = annotation_colors,
        gaps_col = patient_cluster_gaps,
        border_color = NA,
        cellwidth = 0.07,    # Much much smaller cell width
        cellheight = 12,      # Even smaller cell height
        show_rownames = TRUE,
        show_colnames = show_colnames,
        fontsize_row = 12,    # Even smaller font for gene name
        fontsize = 8,        # Smaller general font size
        main = paste0("Pan-Cancer Chaperone Expression - ", k, " Patient Clusters\n",
                     ncol(expression_matrix), " patients from ", 
                     length(combined_data$project_info), " projects"),
        filename = NA,
        na_col = "grey"  # Handle NAs explicitly
      )
    }, error = function(e) {
      message("Error generating heatmap for k=", k, ": ", e$message)
      message("  Expression matrix summary:")
      message("    Dimensions: ", paste(dim(expression_matrix), collapse = " x "))
      message("    Class: ", class(expression_matrix))
      message("    Range: ", paste(range(expression_matrix, na.rm = TRUE), collapse = " to "))
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), cex = 1.2)
    })
    
    dev.off()
    message("Generated heatmap with ", k, " clusters: ", heatmap_file)
    
    # Save patient cluster information for 7 clusters
    if (k == 7) {
      # Create patient ID lists for each cluster
      patient_cluster_lists <- list()
      
      for (cluster_num in 1:k) {
        cluster_patients <- patient_metadata$Patient_ID[patient_clusters == cluster_num]
        patient_cluster_lists[[paste0("Cluster_", cluster_num)]] <- cluster_patients
      }
      
      # Save patient IDs by cluster to separate files for easy manual stratification
      for (cluster_num in 1:k) {
        cluster_file <- file.path(pictures_dir, paste0("cluster_", cluster_num, "_patient_ids.txt"))
        cluster_patients <- patient_metadata$Patient_ID[patient_clusters == cluster_num]
        writeLines(cluster_patients, cluster_file)
        message("Saved ", length(cluster_patients), " patient IDs for Cluster ", cluster_num, " to: ", cluster_file)
      }
      
      # Also create a single summary file with all clusters
      cluster_summary_file <- file.path(pictures_dir, "all_clusters_7_patient_ids.csv")
      patient_cluster_mapping <- data.frame(
        Patient_ID = patient_metadata$Patient_ID,
        Cluster = patient_clusters,
        stringsAsFactors = FALSE
      )
      write.csv(patient_cluster_mapping, cluster_summary_file, row.names = FALSE)
      message("Saved complete patient-to-cluster mapping to: ", cluster_summary_file)
    }
  }
  
  message("All dendrograms and heatmaps saved to: ", pictures_dir)
}

# Main execution function
main <- function(force_recompute_data = FALSE, force_recompute_clustering = FALSE) {
  
  # Step 1: Combine all projects (using existing caching)
  combined_data <- combine_all_projects(force_recompute = force_recompute_data)
  
  # Step 2: Perform clustering
  clustering_results <- perform_clustering(combined_data, 
                                         force_recompute = force_recompute_clustering)
  
  # Step 3: Generate dendrogram and heatmaps with tree cuts
  generate_dendrogram_analysis(combined_data, clustering_results)
  
  # Print summary
  message("\n=== Dendrogram Analysis Complete! ===")
  message("Total patients analyzed: ", ncol(combined_data$combined_expression))
  message("Total projects included: ", length(combined_data$project_info))
  message("Common chaperone genes: ", length(combined_data$common_genes))
  
  message("\nGenerated files in ", pictures_dir, ":")
  message("- patient_dendrogram.pdf")
  message("- heatmap_3_clusters.pdf")
  message("- heatmap_5_clusters.pdf")
  message("- heatmap_7_clusters.pdf")
  message("- cluster_1_patient_ids.txt through cluster_7_patient_ids.txt (patient IDs for each cluster)")
  message("- all_clusters_7_patient_ids.csv (complete patient-to-cluster mapping)")
  
  message("\nCached data available for future runs:")
  message("- Combined expression data: ", combined_data_cache)
  message("- Clustering results: ", clustering_cache)
  
  return(list(combined_data = combined_data, clustering_results = clustering_results))
}

# Run the analysis
# Set force_recompute_data=TRUE to reload and reprocess all projects from scratch
# Set force_recompute_clustering=TRUE to recompute clustering
result <- main(force_recompute_data = FALSE, force_recompute_clustering = FALSE)

message("\n=== Dendrogram analysis completed at: ", Sys.time(), " ===")
