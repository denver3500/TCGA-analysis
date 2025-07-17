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

# Define file paths
project_counts_file <- "TCGA-Chaperones/higher_in_tumor_gene_counts.csv"
gene_list_file <- "TCGA-Chaperones/gene_list.csv"
rds_dir <- "TCGA-Chaperones/rds"
output_dir <- "TCGA-Chaperones/heatmap_patients_stratification/deseq2_clustering"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load project data and filter for projects with >= 10 higher_in_tumor_count
projects_data <- read_csv(project_counts_file, show_col_types = FALSE) %>%
  filter(higher_in_tumor_count >= 10) 

# Load gene list
gene_list <- read_csv(gene_list_file, show_col_types = FALSE)
message("Loaded gene list with ", nrow(gene_list), " genes")

# Function to extract base ENSEMBL ID without version
extract_base_ensembl <- function(ensembl_id) {
  sapply(strsplit(ensembl_id, "\\."), function(x) x[1])
}

# Function to determine optimal number of clusters
determine_optimal_clusters <- function(data, max_k = 10, method = "silhouette") {
  if (method == "silhouette") {
    # Silhouette method
    sil_width <- c()
    for (i in 2:max_k) {
      km <- kmeans(data, centers = i, nstart = 25)
      ss <- silhouette(km$cluster, dist(data))
      sil_width[i-1] <- mean(ss[,3])
    }
    k_optimal <- which.max(sil_width) + 1
  } else {
    # Default: use fixed number based on dataset size
    n_samples <- ncol(data)
    if (n_samples < 30) {
      k_optimal <- 2
    } else if (n_samples < 100) {
      k_optimal <- 3 
    } else if (n_samples < 200) {
      k_optimal <- 4
    } else {
      k_optimal <- 5
    }
  }
  return(k_optimal)
}

# Main function to process each project following DESeq2 best practices
process_project_deseq2 <- function(project_id) {
  message("\n=== Processing project with DESeq2: ", project_id, " ===")
  
  # Find and load RDS file
  rds_files <- list.files(rds_dir, pattern = paste0("^", project_id, "_.*\\.rds$"), full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("  No RDS file found for project ", project_id)
    return(NULL)
  }
  
  message("  Loading RDS file: ", basename(rds_files[1]))
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
  message("  Using ", n_samples, " tumor samples")
  
  # Get RAW COUNT data for DESeq2 (required!)
  assay_names <- assayNames(se_object)
  
  # DESeq2 needs raw counts
  if ("unstranded" %in% assay_names) {
    count_data <- assay(se_object, "unstranded")
    used_assay <- "unstranded"
    message("  Using unstranded count data for DESeq2")
  } else if ("stranded_first" %in% assay_names) {
    count_data <- assay(se_object, "stranded_first")
    used_assay <- "stranded_first"
    message("  Using stranded_first count data for DESeq2")
  } else if ("stranded_second" %in% assay_names) {
    count_data <- assay(se_object, "stranded_second")
    used_assay <- "stranded_second"
    message("  Using stranded_second count data for DESeq2")
  } else {
    # If no raw counts available, fall back to first assay
    count_data <- assay(se_object, assay_names[1])
    used_assay <- assay_names[1]
    message("  WARNING: Using ", assay_names[1], " for DESeq2 (may not be raw counts)")
  }
  
  # Round counts to integers (DESeq2 requirement)
  count_data <- round(count_data)
  
  message("  Full count matrix: ", nrow(count_data), " genes x ", ncol(count_data), " samples")
  
  # Create sample metadata for DESeq2
  sample_metadata <- data.frame(
    sample_id = colnames(count_data),
    sample_type = "tumor",
    stringsAsFactors = FALSE
  )
  rownames(sample_metadata) <- colnames(count_data)
  
  message("  Created sample metadata for ", nrow(sample_metadata), " samples")
  
  # STEP 1: Create DESeq2 dataset with ALL GENES
  message("  Creating DESeq2 dataset with ALL genes...")
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_metadata,
    design = ~ 1  # Intercept only since no comparison groups
  )
  
  message("  DESeq2 dataset created with ", nrow(dds), " genes and ", ncol(dds), " samples")
  
  # STEP 2: Filter out genes with very low counts across the entire dataset
  # Keep genes with at least 10 counts in at least 25% of samples
  min_samples <- ceiling(0.25 * ncol(dds))
  keep <- rowSums(counts(dds) >= 10) >= min_samples
  dds <- dds[keep, ]
  
  message("  After filtering low-count genes: ", nrow(dds), " genes retained")
  
  # STEP 3: Perform DESeq2 normalization on ALL GENES
  message("  Running DESeq2 normalization on all genes...")
  dds <- DESeq(dds)
  
  # STEP 4: Apply VST transformation to ALL GENES
  message("  Applying variance stabilizing transformation to all genes...")
  vsd <- vst(dds, blind = TRUE)
  vst_all_genes <- assay(vsd)
  
  message("  VST transformation complete!")
  message("    VST data range: ", round(min(vst_all_genes), 2), " to ", round(max(vst_all_genes), 2))
  
  # STEP 5: Extract chaperone genes from transformed data
  message("\n  === Extracting chaperone genes from transformed data ===")
  
  # Get all gene IDs and create mapping
  all_genes <- rownames(vst_all_genes)
  
  if (any(grepl("\\.", all_genes))) {
    # Remove version numbers from ENSEMBL IDs
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
  
  # Find chaperone genes in the transformed dataset
  chaperone_genes <- c()
  chaperone_names <- c()
  chaperone_families <- c()
  
  for (i in 1:nrow(gene_list)) {
    ensembl_id <- gene_list$ENSEMBL[i]
    gene_name <- gene_list$Name[i]
    family <- gene_list$Family[i]
    
    # Find matching gene in dataset
    match_idx <- which(gene_mapping$base_id == ensembl_id)
    
    if (length(match_idx) > 0) {
      full_id <- gene_mapping$full_id[match_idx[1]]
      chaperone_genes <- c(chaperone_genes, full_id)
      chaperone_names <- c(chaperone_names, gene_name)
      chaperone_families <- c(chaperone_families, family)
    }
  }
  
  if (length(chaperone_genes) == 0) {
    message("  No chaperone genes found in transformed dataset")
    return(NULL)
  }
  
  message("  Found ", length(chaperone_genes), " of ", nrow(gene_list), " chaperone genes in transformed data")
  
  # Extract chaperone expression data
  vst_chaperones <- vst_all_genes[chaperone_genes, , drop = FALSE]
  rownames(vst_chaperones) <- chaperone_names
  
  message("  Chaperone expression matrix: ", nrow(vst_chaperones), " genes x ", ncol(vst_chaperones), " patients")
  
  # STEP 6: BIDIRECTIONAL CLUSTERING ANALYSIS
  message("\n  === Performing bidirectional clustering ===")
  
  # Set up output subdirectory for this project
  project_dir <- file.path(output_dir, project_id)
  dir.create(project_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Compute distance matrices for genes and patients
  gene_dist <- dist(vst_chaperones)
  patient_dist <- dist(t(vst_chaperones))
  
  # Perform hierarchical clustering
  gene_hclust <- hclust(gene_dist, method = "ward.D2")
  patient_hclust <- hclust(patient_dist, method = "ward.D2")
  
  # Create dendrograms FIRST before using them
  gene_dend <- as.dendrogram(gene_hclust)
  patient_dend <- as.dendrogram(patient_hclust)
  
  # Determine optimal number of clusters based on dataset size
  k_genes <- min(nrow(vst_chaperones), 5)  # Gene clusters (max 5)
  
  # Try to estimate optimal number of patient clusters
  if(n_samples >= 20) {
    # Use silhouette score if enough samples
    k_patients <- determine_optimal_clusters(t(vst_chaperones), max_k = min(8, floor(n_samples/10)))
    message("  Determined optimal number of patient clusters: ", k_patients)
  } else {
    # For small datasets, use simple rule
    k_patients <- min(max(2, round(n_samples/15)), 5)
    message("  Using ", k_patients, " patient clusters based on sample size")
  }
  
  # Cut dendrograms to get clusters
  gene_clusters <- cutree(gene_hclust, k = k_genes)
  patient_clusters <- cutree(patient_hclust, k = k_patients)
  
  message("  Created ", k_genes, " gene clusters and ", k_patients, " patient clusters")
  
  # Create annotation data frames
  gene_annotation <- data.frame(
    Cluster = factor(gene_clusters),
    Family = factor(chaperone_families)
  )
  rownames(gene_annotation) <- rownames(vst_chaperones)
  
  patient_annotation <- data.frame(
    Cluster = factor(patient_clusters)
  )
  rownames(patient_annotation) <- colnames(vst_chaperones)
  
  # Set up colors for annotations
  gene_colors <- list(
    Cluster = setNames(viridis(k_genes), levels(factor(gene_clusters))),
    Family = setNames(rainbow(length(unique(chaperone_families))), 
                     levels(factor(chaperone_families)))
  )
  
  patient_colors <- list(
    Cluster = setNames(plasma(k_patients), levels(factor(patient_clusters)))
  )
  
  # Save clustering objects for later re-analysis
  clustering_objects <- list(
    gene_dist = gene_dist,
    patient_dist = patient_dist,
    gene_hclust = gene_hclust,
    patient_hclust = patient_hclust,
    gene_dend = gene_dend,       # Now these are defined before being used
    patient_dend = patient_dend, # Now these are defined before being used
    vst_chaperones = vst_chaperones,
    gene_clusters = gene_clusters,
    patient_clusters = patient_clusters,
    k_genes = k_genes,
    k_patients = k_patients
  )

  # Save as RDS file for later manipulation
  clustering_rds_file <- file.path(project_dir, paste0(project_id, "_clustering_objects.rds"))
  saveRDS(clustering_objects, clustering_rds_file)
  message("  Saved clustering objects for post-hoc analysis: ", clustering_rds_file)

  # Generate enhanced heatmap with annotations
  message("  Generating enhanced clustered heatmap...")
  
  enhanced_heatmap_file <- file.path(project_dir, paste0(project_id, "_bidirectional_clustered_heatmap.pdf"))
  pdf(enhanced_heatmap_file, width = 16, height = 10)
  
  pheatmap_result <- pheatmap(
    vst_chaperones,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_rows = gene_hclust,
    cluster_cols = patient_hclust,
    annotation_row = gene_annotation,
    annotation_col = patient_annotation,
    annotation_colors = c(gene_colors, patient_colors),
    show_rownames = TRUE,
    show_colnames = ncol(vst_chaperones) <= 60,
    fontsize_row = 8,
    fontsize_col = 6,
    main = paste0(project_id, " - Bidirectional Clustering: ", 
                 k_genes, " gene clusters × ", k_patients, " patient clusters"),
    filename = NA
  )
  
  dev.off()
  message("  Generated enhanced clustered heatmap: ", enhanced_heatmap_file)
  
  # Save dendrograms separately - no need to recreate them as we already have them
  dend_file <- file.path(project_dir, paste0(project_id, "_dendrograms.pdf"))
  pdf(dend_file, width = 12, height = 10)
  
  # Plot patient dendrogram with colored clusters
  par(mar = c(4, 2, 3, 8))
  patient_dend_colored <- color_branches(patient_dend, k = k_patients)
  plot(patient_dend_colored, main = paste0(project_id, " - Patient Clustering Dendrogram"), 
       horiz = TRUE, axes = FALSE)
  colored_labels <- cutree(patient_hclust, k = k_patients)
  legend("topleft", legend = paste("Cluster", 1:k_patients), 
         fill = plasma(k_patients), title = "Patient Clusters", cex = 0.8)
  
  # Plot gene dendrogram with colored clusters
  par(mar = c(4, 2, 3, 8))
  gene_dend_colored <- color_branches(gene_dend, k = k_genes)
  plot(gene_dend_colored, main = paste0(project_id, " - Gene Clustering Dendrogram"),  
       horiz = TRUE, axes = FALSE)
  legend("topleft", legend = paste("Cluster", 1:k_genes), 
         fill = viridis(k_genes), title = "Gene Clusters", cex = 0.8)
  
  dev.off()
  message("  Generated dendrograms: ", dend_file)
  
  # Save cluster assignments
  cluster_assignments <- data.frame(
    Patient_ID = names(patient_clusters),
    Patient_Cluster = patient_clusters,
    stringsAsFactors = FALSE
  )
  
  cluster_file <- file.path(project_dir, paste0(project_id, "_patient_clusters.csv"))
  write_csv(cluster_assignments, cluster_file)
  message("  Saved patient cluster assignments: ", cluster_file)
  
  gene_clusters_df <- data.frame(
    Gene = names(gene_clusters),
    Gene_Cluster = gene_clusters,
    Family = chaperone_families,
    stringsAsFactors = FALSE
  )
  
  gene_cluster_file <- file.path(project_dir, paste0(project_id, "_gene_clusters.csv"))
  write_csv(gene_clusters_df, gene_cluster_file)
  message("  Saved gene cluster assignments: ", gene_cluster_file)
  
  # Calculate cluster statistics
  cluster_stats <- data.frame()
  
  for(i in 1:k_patients) {
    cluster_patients <- names(patient_clusters[patient_clusters == i])
    cluster_expr <- vst_chaperones[, cluster_patients, drop = FALSE]
    
    cluster_stats <- rbind(cluster_stats, data.frame(
      Project = project_id,
      Cluster = i,
      n_patients = length(cluster_patients),
      pct_of_total = round(100 * length(cluster_patients) / length(patient_clusters), 1),
      mean_expr = round(mean(cluster_expr), 3),
      median_expr = round(median(cluster_expr), 3),
      min_expr = round(min(cluster_expr), 3),
      max_expr = round(max(cluster_expr), 3),
      sd_expr = round(sd(cluster_expr), 3)
    ))
  }
  
  cluster_stats_file <- file.path(project_dir, paste0(project_id, "_cluster_statistics.csv"))
  write_csv(cluster_stats, cluster_stats_file)
  message("  Saved cluster statistics: ", cluster_stats_file)
  
  # Save expression data with cluster assignments
  expression_with_clusters <- vst_chaperones
  colnames(expression_with_clusters) <- paste0(colnames(vst_chaperones), 
                                             "_C", patient_clusters[colnames(vst_chaperones)])
  
  expr_cluster_file <- file.path(project_dir, paste0(project_id, "_expression_with_clusters.csv"))
  write_csv(as.data.frame(expression_with_clusters) %>% 
             rownames_to_column("Gene"), expr_cluster_file)
  
  message("  Saved expression data with cluster assignments: ", expr_cluster_file)
  
  # Create project summary
  summary_stats <- data.frame(
    Project = project_id,
    n_patients = ncol(vst_chaperones),
    n_chaperone_genes = nrow(vst_chaperones),
    n_total_genes_processed = nrow(dds),
    n_patient_clusters = k_patients,
    n_gene_clusters = k_genes,
    assay_used = used_assay,
    stringsAsFactors = FALSE
  )
  
  # Save summary
  summary_file <- file.path(project_dir, paste0(project_id, "_clustering_summary.csv"))
  write_csv(summary_stats, summary_file)
  message("  Saved project summary: ", summary_file)
  
  message("\n  === Clustering Results ===")
  message("    Patient clusters: ", k_patients, " clusters for ", ncol(vst_chaperones), " patients")
  print(table(patient_clusters))
  message("    Gene clusters: ", k_genes, " clusters for ", nrow(vst_chaperones), " genes")
  print(table(gene_clusters))
  
  return(summary_stats)
}

# Process all projects
message("\n=== Starting bidirectional clustering analysis ===")

results <- list()

for (i in 1:nrow(projects_data)) {
  project_id <- projects_data$project_id[i]
  result <- process_project_deseq2(project_id)
  if (!is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Save overall summary
if (length(results) > 0) {
  summary_df <- bind_rows(results)
  summary_file <- file.path(output_dir, "clustering_summary_all_projects.csv")
  write_csv(summary_df, summary_file)
  
  message("\n=== Clustering Analysis Summary ===")
  print(summary_df)
  message("\nSummary saved to: ", summary_file)
} else {
  message("\nNo projects were successfully processed")
}

message("\n=== Bidirectional clustering analysis complete! ===")
message("\nKey outputs for each project:")
message("  1. Bidirectional clustered heatmap with annotations")
message("  2. Separate dendrograms for genes and patients")
message("  3. Patient cluster assignments and statistics")
message("  4. Gene cluster assignments by family")
message("  5. Expression data with cluster labels")