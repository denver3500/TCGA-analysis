# This script creates three types of correlation visualizations:
# 1. Mean correlation across all projects
# 2. Project similarity clustering
# 3. 4x8 grid layout with corrplot

target_dir <- "TCGA-Chaperones/Gene_Correlation"
if (!endsWith(getwd(), target_dir)) {
  setwd(target_dir)
}

library(tidyverse)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# Set up logging
log_file <- "comprehensive_correlation_analysis_log.txt"
con <- file(log_file, "w")
sink(con, split = TRUE)

message("=== COMPREHENSIVE TCGA CORRELATION ANALYSIS ===")
message("Starting at: ", Sys.time())
message("R version: ", R.version.string)
message("Working directory: ", getwd())

# Define order for 4x8 grid
project_order <- c(
  "TCGA-DLBC", "TCGA-PRAD", "TCGA-TGCT", "TCGA-THYM",
  "TCGA-CHOL", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
  "TCGA-LGG", "TCGA-MESO", "TCGA-PCPG", "TCGA-THCA",
  "TCGA-SKCM", "TCGA-UVM", "TCGA-ACC", "TCGA-BLCA",
  "TCGA-BRCA", "TCGA-CESC", "TCGA-COAD", "TCGA-ESCA",
  "TCGA-GBM", "TCGA-HNSC", "TCGA-LIHC", "TCGA-LUAD",
  "TCGA-LUSC", "TCGA-OV", "TCGA-PAAD", "TCGA-READ",
  "TCGA-SARC", "TCGA-STAD", "TCGA-UCEC", "TCGA-UCS"
)

message("Project order defined for 4x8 grid (32 projects):")
for (i in seq_along(project_order)) {
  message(sprintf("  %d. %s", i, project_order[i]))
}

# Load gene list
gene_list <- read.csv("../gene_list.csv", stringsAsFactors = FALSE)
message("Loaded ", nrow(gene_list), " chaperone genes")
message("Gene names: ", paste(gene_list$Name, collapse = ", "))

# Get project information for better labeling
message("\n=== Loading TCGA project information ===")
projects_info <- getGDCprojects()
projects_info <- projects_info %>%
  filter(str_starts(project_id, "TCGA")) %>%
  select(project_id, name)
message("Found ", nrow(projects_info), " TCGA projects in GDC database")

# Get all TCGA .rds files from raw_data directory
raw_data_dir <- "../../raw_data"
message("\n=== Scanning for TCGA data files ===")
message("Scanning directory: ", raw_data_dir)
tcga_files <- list.files(raw_data_dir, pattern = "^TCGA-.*_transcriptomic_exp\\.rds$", full.names = TRUE)
message("Found ", length(tcga_files), " TCGA files")

# Create output directory
output_dir <- "pictures/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
message("\nOutput directory: ", output_dir)

# Function to extract project ID from filename
extract_project_id <- function(filename) {
  project_id <- basename(filename) %>%
    str_replace("_transcriptomic_exp\\.rds$", "")
  return(project_id)
}

# Function to get project name
get_project_name <- function(project_id) {
  project_name <- projects_info$name[projects_info$project_id == project_id]
  if (length(project_name) > 0) {
    short_name <- str_replace(project_name, "^(.*?)\\s*-\\s*TCGA$", "\\1")
    return(short_name)
  } else {
    return(str_replace(project_id, "TCGA-", ""))
  }
}

# Function to process a single TCGA dataset for correlation analysis
process_tcga_correlation <- function(file_path) {
  project_id <- extract_project_id(file_path)
  message("\n=== Processing ", project_id, " ===")
  
  tryCatch({
    # Load the data
    message("  Loading data...")
    brca.exp <- readRDS(file_path)
    
    # Standardize assay names
    if("unstranded" %in% assayNames(brca.exp)) {
      assayNames(brca.exp)[assayNames(brca.exp) == "unstranded"] <- "counts"
    }
    
    # Get sample metadata
    sample_info <- as.data.frame(colData(brca.exp))
    
    # Filter for tumor samples
    sample.idx <- logical(ncol(brca.exp))
    
    if ("definition" %in% colnames(sample_info)) {
      sample.idx <- sample_info$definition == "Primary solid Tumor"
    } else if ("sample_type" %in% colnames(sample_info)) {
      sample.idx <- sample_info$sample_type == "Primary Tumor"
    } else if ("shortLetterCode" %in% colnames(sample_info)) {
      sample.idx <- sample_info$shortLetterCode == "TP"
    } else {
      message("  ERROR: No valid field found for filtering tumor samples")
      return(NULL)
    }
    
    tumor_samples <- sum(sample.idx)
    message("  Tumor samples: ", tumor_samples)
    
    # Check if we have enough tumor samples
    if (tumor_samples < 10) {
      message("  WARNING: Less than 10 tumor samples. Skipping.")
      return(NULL)
    }
    
    # Subset for tumor samples
    brca.exp.subset <- brca.exp[, sample.idx]
    
    # Extract gene metadata
    gene_metadata <- as.matrix(mcols(rowRanges(brca.exp))[, c("gene_id", "gene_name")])
    
    # Get gene IDs for chaperone genes
    gene_names <- gene_list$Name
    gene_ids <- gene_metadata[gene_metadata[, "gene_name"] %in% gene_names, "gene_id"]
    
    if (length(gene_ids) == 0) {
      message("  ERROR: No matching chaperone gene IDs found")
      return(NULL)
    }
    
    # Subset data for selected genes
    brca.exp.genes <- brca.exp.subset[gene_ids, ]
    
    # Check if we have enough genes
    if (nrow(brca.exp.genes) < 5) {
      message("  ERROR: Less than 5 chaperone genes found")
      return(NULL)
    }
    
    # Create DESeq2 dataset
    message("  Creating DESeq2 dataset...")
    
    # Create simple design matrix
    colData_simple <- data.frame(
      sample = colnames(brca.exp.genes),
      condition = rep("tumor", ncol(brca.exp.genes)),
      row.names = colnames(brca.exp.genes)
    )
    
    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(
      countData = assay(brca.exp.genes, "counts"),
      colData = colData_simple,
      design = ~ 1
    )
    
    # Perform normalization
    message("  Performing DESeq2 normalization...")
    dds <- estimateSizeFactors(dds)
    
    # Get normalized counts and log2 transform
    normalized_counts <- counts(dds, normalized = TRUE)
    log2_counts <- log2(normalized_counts + 1)
    
    # Create gene name mapping
    gene_id_to_name <- setNames(gene_metadata[, "gene_name"], gene_metadata[, "gene_id"])
    
    # Map gene IDs to names
    gene_names_for_plot <- sapply(rownames(log2_counts), function(id) {
      if (id %in% names(gene_id_to_name)) {
        return(gene_id_to_name[id])
      } else {
        return(id)
      }
    })
    
    # Set gene names as rownames
    rownames(log2_counts) <- gene_names_for_plot
    
    # Compute correlation matrix
    message("  Computing correlation matrix...")
    cor_matrix <- cor(t(log2_counts), method = "pearson")
    
    message("  SUCCESS: ", nrow(log2_counts), " genes, ", ncol(log2_counts), " samples")
    
    return(list(
      project_id = project_id,
      correlation_matrix = cor_matrix,
      gene_names = gene_names_for_plot,
      n_samples = ncol(log2_counts),
      n_genes = nrow(log2_counts)
    ))
    
  }, error = function(e) {
    message("  ERROR: Failed to process ", project_id, ": ", e$message)
    return(NULL)
  })
}

# Process all TCGA files
message("\n=== PROCESSING ALL TCGA DATASETS ===")
correlation_results <- map(tcga_files, process_tcga_correlation)
names(correlation_results) <- map_chr(tcga_files, extract_project_id)

# Remove NULL results
correlation_results <- correlation_results[!map_lgl(correlation_results, is.null)]
message("Successfully processed ", length(correlation_results), " TCGA datasets")

# Get common genes across all projects
all_genes <- unique(unlist(map(correlation_results, ~ rownames(.x$correlation_matrix))))
correlation_matrices_list <- map(correlation_results, ~ .x$correlation_matrix)
common_genes <- Reduce(intersect, map(correlation_matrices_list, rownames))
message("Common genes across all projects: ", length(common_genes))

# Standardize all correlation matrices
standardized_matrices <- map(correlation_matrices_list, function(mat) {
  mat[common_genes, common_genes]
})

# Create project names for labeling
project_labels <- sapply(names(standardized_matrices), get_project_name)

# ============================================================================
# 1. MEAN CORRELATION ACROSS ALL PROJECTS
# ============================================================================

message("\n=== CREATING MEAN CORRELATION HEATMAP ===")

# Calculate mean correlations across all projects
mean_correlations <- map_dfr(standardized_matrices, function(mat) {
  expand_grid(gene1 = rownames(mat), gene2 = colnames(mat)) %>%
    mutate(correlation = as.vector(mat)) %>%
    filter(gene1 != gene2)  # Remove self-correlations
}) %>%
  group_by(gene1, gene2) %>%
  summarise(
    mean_correlation = mean(correlation, na.rm = TRUE),
    n_projects = n(),
    .groups = "drop"
  )

# Create mean correlation matrix
mean_cor_matrix <- matrix(NA, nrow = length(common_genes), ncol = length(common_genes))
rownames(mean_cor_matrix) <- colnames(mean_cor_matrix) <- common_genes
diag(mean_cor_matrix) <- 1

# Fill the matrix
for(i in 1:nrow(mean_correlations)) {
  gene1 <- mean_correlations$gene1[i]
  gene2 <- mean_correlations$gene2[i]
  corr_val <- mean_correlations$mean_correlation[i]
  
  mean_cor_matrix[gene1, gene2] <- corr_val
  mean_cor_matrix[gene2, gene1] <- corr_val
}

# Create mean correlation heatmap
png(file.path(output_dir, "mean_correlation_heatmap.png"), 
    width = 14, height = 12, units = "in", res = 300, bg = "white")

# Color
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Heatmap
ht <- Heatmap(mean_cor_matrix,
              col = col_fun,
              name = "Correlation",  # This sets the legend title
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 12),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(!is.na(mean_cor_matrix[i, j]))
                  grid.text(sprintf("%.2f", mean_cor_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              column_title = "Mean mitochondrial chaperone correlations across different TCGA projects",
              column_title_gp = gpar(fontsize = 16, fontface = "bold"),
              heatmap_legend_param = list(
                title = "Correlation",
                title_gp = gpar(fontsize = 14, fontface = "bold"),
                labels_gp = gpar(fontsize = 12)
              ))

draw(ht)

dev.off()
message("Saved mean correlation heatmap: mean_correlation_heatmap.png")

# ============================================================================
# 2. PROJECT SIMILARITY CLUSTERING
# ============================================================================

message("\n=== CREATING PROJECT SIMILARITY CLUSTERING ===")

# Calculate distance matrix for project clustering
project_distances <- matrix(0, nrow = length(standardized_matrices), ncol = length(standardized_matrices))
rownames(project_distances) <- colnames(project_distances) <- names(standardized_matrices)

for (i in 1:length(standardized_matrices)) {
  for (j in 1:length(standardized_matrices)) {
    if (i != j) {
      mat1 <- standardized_matrices[[i]]
      mat2 <- standardized_matrices[[j]]
      
      # Get upper triangle correlations
      upper_tri <- upper.tri(mat1)
      corr1 <- mat1[upper_tri]
      corr2 <- mat2[upper_tri]
      
      # Calculate Euclidean distance between correlation patterns
      project_distances[i, j] <- sqrt(sum((corr1 - corr2)^2))
    }
  }
}

# Perform hierarchical clustering
project_hclust <- hclust(as.dist(project_distances), method = "ward.D2")

# Create similarity matrix (inverse of distance)
max_dist <- max(project_distances)
project_similarity <- max_dist - project_distances

# Create clustered similarity heatmap
png(file.path(output_dir, "project_clustering.png"), 
    width = 20, height = 16, units = "in", res = 300, bg = "white")

# Get project names for labeling
project_names_full <- sapply(rownames(project_similarity), get_project_name)

# Color
sim_col_fun <- colorRamp2(c(min(project_similarity), 
                           min(project_similarity) + (max(project_similarity) - min(project_similarity)) * 1/3,
                           min(project_similarity) + (max(project_similarity) - min(project_similarity)) * 2/3,
                           max(project_similarity)), 
                         c("#fff0f3", "#ff8fa3", "#c9184a", "#590d22"))

# Heatmap
ht2 <- Heatmap(project_similarity,
               col = sim_col_fun,
               name = "Similarity",  # This sets the legend title
               cluster_rows = project_hclust,
               cluster_columns = project_hclust,
               show_row_names = TRUE,
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 12),
               column_names_gp = gpar(fontsize = 12),
               row_labels = project_names_full,
               column_labels = project_names_full,
               column_title = "Mitochondrial chaperone correlation similarity in different TCGA projects",
               column_title_gp = gpar(fontsize = 14, fontface = "bold"),
               # Move row names to left and dendrogram to right
               row_names_side = "left",
               row_dend_side = "right",
               # Rotate column names to 45 degrees
               column_names_rot = 45,
               # Adjust margins to prevent overlap
               left_annotation = NULL,
               right_annotation = NULL,
               # Add margins for long labels
               width = unit(10, "in"),
               height = unit(10, "in"),
               heatmap_legend_param = list(
                 title = "Similarity",
                 title_gp = gpar(fontsize = 14, fontface = "bold"),
                 labels_gp = gpar(fontsize = 12),
                 # Position legend to avoid overlap
                 legend_direction = "vertical",
                 legend_width = unit(4, "cm")
               ))

# Draw with padding to ensure labels fit
draw(ht2, padding = unit(c(3, 3, 3, 3), "cm"))

dev.off()
message("Saved project similarity clustering: project_similarity_clustering.png")

# ============================================================================
# 3. 4x8 GRID LAYOUT WITH CORRPLOT
# ============================================================================

message("\n=== CREATING 4x8 GRID CORRELATION PLOTS ===")

# Filter for projects in our predefined order
available_projects <- intersect(project_order, names(standardized_matrices))
message("Projects available for 4x8 grid: ", length(available_projects))

# Create ordered list of matrices
ordered_matrices <- list()
ordered_labels <- character()

for (project_id in project_order) {
  if (project_id %in% names(standardized_matrices)) {
    ordered_matrices[[project_id]] <- standardized_matrices[[project_id]]
    ordered_labels[project_id] <- project_labels[project_id]
  }
}

# Sort genes alphabetically
genes_sorted <- sort(common_genes)

# Create corrplot version
message("Creating corrplot 4x8 grid...")

png(file.path(output_dir, "tcga_correlations_4x8_corrplot.png"), 
    width = 20, height = 26, units = "in", res = 300, bg = "white")

# Set up 4x8 grid layout
par(mfrow = c(8, 4), mar = c(3, 3, 3, 1), oma = c(6, 6, 8, 6))

# Plot matrices in the specified order
for (i in 1:length(ordered_matrices)) {
  project_id <- names(ordered_matrices)[i]
  project_name <- ordered_labels[project_id]
  mat <- ordered_matrices[[project_id]][genes_sorted, genes_sorted]
  
  corrplot(mat, 
           method = "color",
           type = "full", 
           order = "original",
           col = colorRampPalette(c("blue", "white", "red"))(100),
           tl.cex = 0.8,
           tl.col = "black", 
           tl.srt = 45,
           cl.pos = "r",  
           cl.cex = 0.6,
           title = project_name,
           mar = c(1, 1, 2, 1))
}

# Add main title
mtext("Mitochondrial chaperone gene correlations across TCGA projects", 
      side = 3, outer = TRUE, cex = 1.8, font = 2, line = 4)

dev.off()
message("Saved 4x8 corrplot grid: tcga_correlations_4x8_corrplot.png")

# ============================================================================
# SAVE RESULTS AND SUMMARY
# ============================================================================

message("\n=== SAVING RESULTS ===")

# Save all correlation data
saveRDS(list(
  ordered_matrices = ordered_matrices,
  project_labels = ordered_labels,
  common_genes = common_genes,
  mean_correlation_matrix = mean_cor_matrix,
  project_similarity_matrix = project_similarity,
  project_clustering = project_hclust
), file.path(output_dir, "comprehensive_correlation_results.rds"))

# Save mean correlations
write_csv(mean_correlations, file.path(output_dir, "mean_correlations.csv"))

# Save project summary
project_summary <- data.frame(
  project_id = names(correlation_results),
  project_name = sapply(names(correlation_results), get_project_name),
  n_genes = map_int(correlation_results, ~ .x$n_genes),
  n_samples = map_int(correlation_results, ~ .x$n_samples),
  included_in_grid = names(correlation_results) %in% available_projects,
  stringsAsFactors = FALSE
)

write_csv(project_summary, file.path(output_dir, "project_summary.csv"))

# ============================================================================
# FINAL SUMMARY
# ============================================================================

message("\n=== ANALYSIS COMPLETE ===")
message("Successfully created three types of correlation visualizations:")
message("1. Mean correlation heatmap: mean_correlation_heatmap.png")
message("2. Project similarity clustering: project_similarity_clustering.png")
message("3. 4x8 grid corrplot: tcga_correlations_4x8_corrplot.png")
message("3. 4x8 grid corrplot: tcga_correlations_4x8_corrplot.png")
message("")
message("Data processed:")
message("  - TCGA projects: ", length(correlation_results))
message("  - Common genes: ", length(common_genes))
message("  - Projects in 4x8 grid: ", length(available_projects))
message("  - Normalization: DESeq2 with log2 transformation")
message("")
message("All files saved to: ", output_dir)
message("Analysis completed at: ", Sys.time())

# Close logging
sink()
close(con)

message("Comprehensive correlation analysis completed successfully!")
