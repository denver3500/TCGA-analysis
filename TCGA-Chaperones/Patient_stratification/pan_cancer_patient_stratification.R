# Pan-Cancer Patient Stratification Based on Chaperone Expression
# This script loads all TCGA projects at once and stratifies all patients together
# based on chaperone expression patterns using DESeq2 normalization

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

# Define file paths - adjust relative to current directory
gene_list_file <- "../gene_list.csv"
raw_data_dir <- "../../raw_data"
output_dir <- "pan_cancer_clustering"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define cache files
cache_dir <- file.path(output_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

combined_data_cache <- file.path(cache_dir, "combined_expression_data.rds")
clustering_cache <- file.path(cache_dir, "clustering_results.rds")

message("=== Pan-Cancer Patient Stratification Based on Chaperone Expression ===")
message("Starting at: ", Sys.time())

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
}

# Function to determine optimal number of clusters
determine_optimal_clusters <- function(data, max_k = 15, method = "silhouette") {
  if (method == "silhouette") {
    message("    Determining optimal clusters using silhouette method...")
    sil_width <- c()
    for (i in 2:max_k) {
      if (i %% 5 == 0) message("      Testing k=", i, "...")
      tryCatch({
        km <- kmeans(data, centers = i, nstart = 25, iter.max = 100)
        ss <- silhouette(km$cluster, dist(data))
        sil_width[i-1] <- mean(ss[,3])
      }, error = function(e) {
        sil_width[i-1] <- 0
      })
    }
    k_optimal <- which.max(sil_width) + 1
    message("    Optimal k found: ", k_optimal, " (silhouette score: ", round(max(sil_width), 3), ")")
  } else {
    # Default: use fixed number based on dataset size
    n_samples <- ncol(data)
    if (n_samples < 100) {
      k_optimal <- 3
    } else if (n_samples < 500) {
      k_optimal <- 5
    } else if (n_samples < 1000) {
      k_optimal <- 7
    } else if (n_samples < 2000) {
      k_optimal <- 10
    } else {
      k_optimal <- 12
    }
    message("    Using k=", k_optimal, " based on dataset size (", n_samples, " samples)")
  }
  return(k_optimal)
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

# Main function to combine all projects
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

# Function to perform pan-cancer clustering
perform_pan_cancer_clustering <- function(combined_data, force_recompute = FALSE) {
  
  # Check if cached clustering exists
  if (!force_recompute && file.exists(clustering_cache)) {
    message("\n=== Loading cached clustering results ===")
    clustering_results <- readRDS(clustering_cache)
    message("Loaded clustering results for ", length(clustering_results$patient_clusters), " patients")
    return(clustering_results)
  }
  
  message("\n=== Performing pan-cancer patient clustering ===")
  
  expression_matrix <- combined_data$combined_expression
  patient_metadata <- combined_data$patient_metadata
  gene_families <- combined_data$gene_families
  
  # Compute distance matrices
  message("Computing distance matrices...")
  gene_dist <- dist(expression_matrix)
  patient_dist <- dist(t(expression_matrix))
  
  # Perform hierarchical clustering
  message("Performing hierarchical clustering...")
  gene_hclust <- hclust(gene_dist, method = "ward.D2")
  patient_hclust <- hclust(patient_dist, method = "ward.D2")
  
  # Create dendrograms
  gene_dend <- as.dendrogram(gene_hclust)
  patient_dend <- as.dendrogram(patient_hclust)
  
  # Determine optimal clusters
  k_genes <- min(nrow(expression_matrix), 6)  # Max 6 gene clusters
  
  # For patients, use silhouette method for large datasets
  n_patients <- ncol(expression_matrix)
  if (n_patients >= 100) {
    k_patients <- determine_optimal_clusters(t(expression_matrix), 
                                           max_k = min(15, floor(n_patients/50)), 
                                           method = "silhouette")
  } else {
    k_patients <- min(max(3, round(n_patients/20)), 8)
    message("Using k=", k_patients, " patient clusters based on dataset size")
  }
  
  # Cut dendrograms
  gene_clusters <- cutree(gene_hclust, k = k_genes)
  patient_clusters <- cutree(patient_hclust, k = k_patients)
  
  message("Created ", k_genes, " gene clusters and ", k_patients, " patient clusters")
  
  # Create detailed patient metadata with clusters
  enhanced_patient_metadata <- patient_metadata %>%
    mutate(
      Patient_Cluster = patient_clusters[Patient_ID],
      Project_Cluster = paste0(Project, "_C", Patient_Cluster)
    )
  
  # Calculate cluster statistics by project
  cluster_stats <- enhanced_patient_metadata %>%
    group_by(Project, Patient_Cluster) %>%
    summarise(
      n_patients = n(),
      .groups = 'drop'
    ) %>%
    group_by(Project) %>%
    mutate(
      pct_of_project = round(100 * n_patients / sum(n_patients), 1)
    ) %>%
    ungroup()
  
  # Calculate overall cluster statistics
  overall_cluster_stats <- enhanced_patient_metadata %>%
    group_by(Patient_Cluster) %>%
    summarise(
      n_patients = n(),
      n_projects = n_distinct(Project),
      projects = paste(unique(Project), collapse = ", "),
      pct_of_total = round(100 * n_patients / nrow(enhanced_patient_metadata), 1),
      .groups = 'drop'
    ) %>%
    arrange(Patient_Cluster)
  
  # Calculate expression statistics per cluster
  expression_stats <- data.frame()
  for(i in 1:k_patients) {
    cluster_patients <- enhanced_patient_metadata$Patient_ID[enhanced_patient_metadata$Patient_Cluster == i]
    cluster_expr <- expression_matrix[, cluster_patients, drop = FALSE]
    
    expression_stats <- rbind(expression_stats, data.frame(
      Patient_Cluster = i,
      n_patients = length(cluster_patients),
      mean_expr = round(mean(cluster_expr), 3),
      median_expr = round(median(cluster_expr), 3),
      sd_expr = round(sd(cluster_expr), 3),
      min_expr = round(min(cluster_expr), 3),
      max_expr = round(max(cluster_expr), 3)
    ))
  }
  
  # Create gene cluster information
  gene_cluster_info <- data.frame(
    Gene = names(gene_clusters),
    Gene_Cluster = gene_clusters,
    Family = gene_families[names(gene_clusters)],
    stringsAsFactors = FALSE
  )
  
  # Combine all clustering results
  clustering_results <- list(
    # Core clustering objects
    gene_hclust = gene_hclust,
    patient_hclust = patient_hclust,
    gene_dend = gene_dend,
    patient_dend = patient_dend,
    gene_clusters = gene_clusters,
    patient_clusters = patient_clusters,
    k_genes = k_genes,
    k_patients = k_patients,
    
    # Metadata and statistics
    enhanced_patient_metadata = enhanced_patient_metadata,
    gene_cluster_info = gene_cluster_info,
    cluster_stats_by_project = cluster_stats,
    overall_cluster_stats = overall_cluster_stats,
    expression_stats = expression_stats,
    
    # Distance matrices
    gene_dist = gene_dist,
    patient_dist = patient_dist
  )
  
  # Cache clustering results
  message("Caching clustering results to: ", clustering_cache)
  saveRDS(clustering_results, clustering_cache)
  
  return(clustering_results)
}

# Function to generate comprehensive visualizations
generate_visualizations <- function(combined_data, clustering_results) {
  message("\n=== Generating pan-cancer clustering visualizations ===")
  
  expression_matrix <- combined_data$combined_expression
  patient_metadata <- clustering_results$enhanced_patient_metadata
  gene_families <- combined_data$gene_families
  
  # Create annotation data frames
  gene_annotation <- data.frame(
    Cluster = factor(clustering_results$gene_clusters),
    Family = factor(gene_families[rownames(expression_matrix)])
  )
  rownames(gene_annotation) <- rownames(expression_matrix)
  
  # Create comprehensive patient annotation
  patient_annotation <- data.frame(
    Project = factor(patient_metadata$Project),
    Cluster = factor(patient_metadata$Patient_Cluster)
  )
  rownames(patient_annotation) <- patient_metadata$Patient_ID
  
  # Set up colors
  n_projects <- length(unique(patient_metadata$Project))
  n_patient_clusters <- clustering_results$k_patients
  n_gene_clusters <- clustering_results$k_genes
  n_families <- length(unique(gene_families))
  
  colors <- list(
    Project = setNames(rainbow(n_projects), levels(patient_annotation$Project)),
    Cluster = setNames(viridis(n_patient_clusters), levels(patient_annotation$Cluster)),
    Family = setNames(brewer.pal(min(n_families, 11), "Set3")[1:n_families], 
                     levels(gene_annotation$Family))
  )
  
  # Add gene cluster colors
  colors$Cluster_genes <- setNames(plasma(n_gene_clusters), 
                                  levels(gene_annotation$Cluster))
  
  # Update colors list for annotations
  annotation_colors <- list(
    Project = colors$Project,
    Cluster = colors$Cluster,
    Family = colors$Family
  )
  names(annotation_colors)[names(annotation_colors) == "Cluster"] <- "Cluster"
  
  # Generate main pan-cancer heatmap
  message("Generating pan-cancer heatmap...")
  
  main_heatmap_file <- file.path(output_dir, "pan_cancer_chaperone_clustering_heatmap.pdf")
  pdf(main_heatmap_file, width = 20, height = 12)
  
  # Show column names only if reasonable number of patients
  show_colnames <- ncol(expression_matrix) <= 200
  
  pheatmap_result <- pheatmap(
    expression_matrix,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_rows = clustering_results$gene_hclust,
    cluster_cols = clustering_results$patient_hclust,
    annotation_row = gene_annotation,
    annotation_col = patient_annotation,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = show_colnames,
    fontsize_row = 8,
    fontsize_col = 4,
    main = paste0("Pan-Cancer Chaperone Expression Clustering\n",
                 ncol(expression_matrix), " patients from ", 
                 length(combined_data$project_info), " projects | ",
                 clustering_results$k_patients, " patient clusters × ",
                 clustering_results$k_genes, " gene clusters"),
    filename = NA
  )
  
  dev.off()
  message("Generated main heatmap: ", main_heatmap_file)
  
  # Generate dendrograms
  message("Generating dendrograms...")
  
  dend_file <- file.path(output_dir, "pan_cancer_dendrograms.pdf")
  pdf(dend_file, width = 16, height = 12)
  
  # Patient dendrogram
  par(mar = c(6, 4, 4, 8))
  patient_dend_colored <- color_branches(clustering_results$patient_dend, 
                                       k = clustering_results$k_patients)
  plot(patient_dend_colored, 
       main = paste0("Pan-Cancer Patient Clustering Dendrogram (", 
                    ncol(expression_matrix), " patients)"), 
       horiz = TRUE, axes = FALSE)
  legend("topleft", 
         legend = paste("Cluster", 1:clustering_results$k_patients), 
         fill = viridis(clustering_results$k_patients), 
         title = "Patient Clusters", cex = 0.8)
  
  # Gene dendrogram
  par(mar = c(6, 4, 4, 8))
  gene_dend_colored <- color_branches(clustering_results$gene_dend, 
                                    k = clustering_results$k_genes)
  plot(gene_dend_colored, 
       main = paste0("Chaperone Gene Clustering Dendrogram (", 
                    nrow(expression_matrix), " genes)"),  
       horiz = TRUE, axes = FALSE)
  legend("topleft", 
         legend = paste("Cluster", 1:clustering_results$k_genes), 
         fill = plasma(clustering_results$k_genes), 
         title = "Gene Clusters", cex = 0.8)
  
  dev.off()
  message("Generated dendrograms: ", dend_file)
  
  message("All visualizations saved to: ", output_dir)
}

# Function to save all results to CSV
save_results_to_csv <- function(combined_data, clustering_results) {
  message("\n=== Saving results to CSV files ===")
  
  # 1. Patient metadata with clusters
  patient_file <- file.path(output_dir, "pan_cancer_patient_clusters.csv")
  write_csv(clustering_results$enhanced_patient_metadata, patient_file)
  message("Saved patient cluster assignments: ", patient_file)
  
  # 2. Gene cluster assignments
  gene_file <- file.path(output_dir, "pan_cancer_gene_clusters.csv")
  write_csv(clustering_results$gene_cluster_info, gene_file)
  message("Saved gene cluster assignments: ", gene_file)
  
  # 3. Cluster statistics by project
  project_stats_file <- file.path(output_dir, "cluster_statistics_by_project.csv")
  write_csv(clustering_results$cluster_stats_by_project, project_stats_file)
  message("Saved cluster statistics by project: ", project_stats_file)
  
  # 4. Overall cluster statistics
  overall_stats_file <- file.path(output_dir, "overall_cluster_statistics.csv")
  write_csv(clustering_results$overall_cluster_stats, overall_stats_file)
  message("Saved overall cluster statistics: ", overall_stats_file)
  
  # 5. Expression statistics per cluster
  expr_stats_file <- file.path(output_dir, "expression_statistics_per_cluster.csv")
  write_csv(clustering_results$expression_stats, expr_stats_file)
  message("Saved expression statistics per cluster: ", expr_stats_file)
  
  # 6. Combined expression data with cluster labels
  expression_with_clusters <- combined_data$combined_expression
  colnames(expression_with_clusters) <- paste0(
    colnames(expression_with_clusters), 
    "_C", clustering_results$patient_clusters[colnames(expression_with_clusters)]
  )
  
  expr_file <- file.path(output_dir, "pan_cancer_expression_with_clusters.csv")
  write_csv(as.data.frame(expression_with_clusters) %>% 
           rownames_to_column("Gene"), expr_file)
  message("Saved expression data with clusters: ", expr_file)
  
  # 7. Project summary
  project_summary <- do.call(rbind, lapply(combined_data$project_info, function(x) {
    data.frame(
      Project = x$project,
      N_Samples = x$n_samples,
      Used_Assay = x$used_assay,
      stringsAsFactors = FALSE
    )
  }))
  
  project_summary_file <- file.path(output_dir, "project_summary.csv")
  write_csv(project_summary, project_summary_file)
  message("Saved project summary: ", project_summary_file)
}

# Main execution
main <- function(force_recompute_data = FALSE, force_recompute_clustering = FALSE) {
  
  # Step 1: Combine all projects
  combined_data <- combine_all_projects(force_recompute = force_recompute_data)
  
  # Step 2: Perform clustering
  clustering_results <- perform_pan_cancer_clustering(combined_data, 
                                                    force_recompute = force_recompute_clustering)
  
  # Step 3: Generate visualizations
  generate_visualizations(combined_data, clustering_results)
  
  # Step 4: Save results to CSV
  save_results_to_csv(combined_data, clustering_results)
  
  # Print summary
  message("\n=== Pan-Cancer Clustering Analysis Complete! ===")
  message("Total patients analyzed: ", ncol(combined_data$combined_expression))
  message("Total projects included: ", length(combined_data$project_info))
  message("Common chaperone genes: ", length(combined_data$common_genes))
  message("Patient clusters identified: ", clustering_results$k_patients)
  message("Gene clusters identified: ", clustering_results$k_genes)
  
  message("\nCluster distribution across projects:")
  print(clustering_results$overall_cluster_stats)
  
  message("\nResults saved to: ", output_dir)
  message("- Main heatmap: pan_cancer_chaperone_clustering_heatmap.pdf")
  message("- Dendrograms: pan_cancer_dendrograms.pdf")
  message("- Patient clusters: pan_cancer_patient_clusters.csv")
  message("- Gene clusters: pan_cancer_gene_clusters.csv")
  message("- Statistics: cluster_statistics_by_project.csv")
  message("- Expression data: pan_cancer_expression_with_clusters.csv")
  
  message("\nCached data available for future runs:")
  message("- Combined expression data: ", combined_data_cache)
  message("- Clustering results: ", clustering_cache)
  message("\nTo force recomputation, set force_recompute_data=TRUE or force_recompute_clustering=TRUE")
  
  return(list(combined_data = combined_data, clustering_results = clustering_results))
}

# Run the analysis
# Set force_recompute_data=TRUE to reload and reprocess all projects from scratch
# Set force_recompute_clustering=TRUE to recompute clustering with different parameters
result <- main(force_recompute_data = FALSE, force_recompute_clustering = FALSE)

message("\n=== Analysis completed at: ", Sys.time(), " ===")
