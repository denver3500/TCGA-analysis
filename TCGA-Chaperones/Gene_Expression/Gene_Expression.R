# This script loads TCGA .rds files from raw_data directory and performs chaperone gene analysis

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(patchwork)
library(SummarizedExperiment)
library(stringr)
library(dplyr)
library(readr)
library(TCGAbiolinks)
library(DESeq2)

# Set up logging 
log_file <- "TCGA-Chaperones/Gene_Expression/analysis_log.txt"
# Create directory if it doesn't exist
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
con <- file(log_file, "w")
sink(con, split = TRUE)

message("=== TCGA Raw Data Analysis ===")
message("Starting at: ", Sys.time())

# Load gene list
gene_list <- read_csv("TCGA-Chaperones/gene_list.csv", show_col_types = FALSE)
chaperone_genes <- gene_list$Name
message("Loaded ", length(chaperone_genes), " chaperone genes")

# Get project information for better labeling
projects_info <- getGDCprojects()
projects_info <- projects_info %>%
  filter(str_starts(project_id, "TCGA")) %>%
  select(project_id, name)

# Get all TCGA .rds files from raw_data directory
raw_data_dir <- "raw_data"
tcga_files <- list.files(raw_data_dir, pattern = "^TCGA-.*_transcriptomic_exp\\.rds$", full.names = TRUE)
message("Found ", length(tcga_files), " TCGA files")

# Function to extract project ID from filename
extract_project_id <- function(filename) {
  basename(filename) %>%
    str_replace("_transcriptomic_exp\\.rds$", "")
}

# Function to process a single TCGA dataset
process_tcga_dataset <- function(file_path) {
  project_id <- extract_project_id(file_path)
  message("Processing ", project_id, "...")
  
  tryCatch({
    # Load the data
    data <- readRDS(file_path)
    
    # Check if it's a SummarizedExperiment object
    if (inherits(data, "SummarizedExperiment")) {
      # Extract expression matrix and sample info
      expr_matrix <- assay(data)
      sample_info <- colData(data)
      gene_info <- rowData(data)
      
      # Get sample types
      if ("sample_type" %in% colnames(sample_info)) {
        normal_samples <- colnames(expr_matrix)[sample_info$sample_type == "Solid Tissue Normal"]
        tumor_samples <- colnames(expr_matrix)[sample_info$sample_type == "Primary Tumor"]
      } else {
        # Fallback to shortLetterCode if available
        if ("shortLetterCode" %in% colnames(sample_info)) {
          normal_samples <- colnames(expr_matrix)[sample_info$shortLetterCode == "NT"]
          tumor_samples <- colnames(expr_matrix)[sample_info$shortLetterCode == "TP"]
        } else {
          normal_samples <- character(0)
          tumor_samples <- colnames(expr_matrix)
        }
      }
      
      # Filter for chaperone genes
      gene_matches <- character(0)
      
      # Try matching with base ENSEMBL IDs
      ensembl_ids <- str_replace(rownames(expr_matrix), "\\.\\d+$", "")
      gene_matches <- rownames(expr_matrix)[ensembl_ids %in% gene_list$ENSEMBL]
      
      # If no matches, try gene symbols
      if (length(gene_matches) == 0 && "gene_name" %in% colnames(gene_info)) {
        gene_matches <- rownames(expr_matrix)[gene_info$gene_name %in% gene_list$Name]
      }
      
      if (length(gene_matches) == 0) {
        message("  No chaperone genes found in ", project_id, ". Skipping.")
        return(NULL)
      }
      
      message("  Found ", length(gene_matches), " chaperone genes in ", project_id)
      message("  Normal samples: ", length(normal_samples))
      message("  Tumor samples: ", length(tumor_samples))
      
      # Filter expression data for chaperone genes
      chaperone_expr <- expr_matrix[gene_matches, , drop = FALSE]
      
      # Create mapping from ENSEMBL to gene names for labeling
      gene_names <- character(length(gene_matches))
      names(gene_names) <- gene_matches
      
      for (i in seq_along(gene_matches)) {
        gene_id <- gene_matches[i]
        base_ensembl <- str_replace(gene_id, "\\.\\d+$", "")
        
        # Find corresponding gene name
        if ("gene_name" %in% colnames(gene_info)) {
          gene_name <- gene_info$gene_name[rownames(expr_matrix) == gene_id]
          if (length(gene_name) > 0 && !is.na(gene_name[1])) {
            gene_names[gene_id] <- gene_name[1]
          }
        }
        
        # If no gene name found, use gene list
        if (gene_names[gene_id] == "" || is.na(gene_names[gene_id])) {
          name_match <- gene_list$Name[gene_list$ENSEMBL == base_ensembl]
          if (length(name_match) > 0) {
            gene_names[gene_id] <- name_match[1]
          } else {
            gene_names[gene_id] <- base_ensembl
          }
        }
      }
      
      # Create result list
      result <- list(
        project_id = project_id,
        expression_data = chaperone_expr,
        normal_samples = normal_samples,
        tumor_samples = tumor_samples,
        gene_matches = gene_matches,
        gene_names = gene_names,
        total_samples = ncol(chaperone_expr)
      )
      
      return(result)
      
    } else {
      message("  Unknown data structure for ", project_id, ". Skipping.")
      return(NULL)
    }
    
  }, error = function(e) {
    message("  Error processing ", project_id, ": ", e$message)
    return(NULL)
  })
}

# Process all TCGA files
message("\\n=== Processing TCGA datasets ===")
tcga_results <- map(tcga_files, process_tcga_dataset)
names(tcga_results) <- map_chr(tcga_files, extract_project_id)

# Remove NULL results
tcga_results <- tcga_results[!map_lgl(tcga_results, is.null)]

message("Successfully processed ", length(tcga_results), " TCGA datasets")

# Create output directory
output_dir <- "TCGA-Chaperones/Gene_Expression/pictures/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to perform differential expression analysis
perform_de_analysis <- function(data_list) {
  results <- tibble()
  normalized_data <- list()
  
  for (project_id in names(data_list)) {
    dataset <- data_list[[project_id]]
    
    message("  Processing ", project_id, " with DESeq2 normalization...")
    
    tryCatch({
      # Prepare count matrix (ensure it's integers)
      count_matrix <- round(dataset$expression_data)
      
      # Remove genes with all zeros
      keep_genes <- rowSums(count_matrix) > 0
      count_matrix <- count_matrix[keep_genes, , drop = FALSE]
      
      # Create sample metadata
      sample_metadata <- data.frame(
        sample_id = colnames(count_matrix),
        condition = ifelse(colnames(count_matrix) %in% dataset$normal_samples, "normal", "tumor"),
        stringsAsFactors = FALSE
      )
      rownames(sample_metadata) <- sample_metadata$sample_id
      
      # Check if we have both normal and tumor samples
      has_normal <- length(dataset$normal_samples) >= 3
      has_tumor <- length(dataset$tumor_samples) >= 3
      
      if (has_normal && has_tumor) {
        # DESeq2 analysis with normal vs tumor comparison
        message("    Performing DESeq2 analysis for ", project_id, " (normal vs tumor)")
        
        # Create DESeq2 object
        dds <- DESeqDataSetFromMatrix(
          countData = count_matrix,
          colData = sample_metadata,
          design = ~ condition
        )
        
        # Filter low count genes
        keep <- rowSums(counts(dds) >= 10) >= 3
        dds <- dds[keep, ]
        
        # Run DESeq2
        dds <- DESeq(dds)
        
        # Get normalized counts
        normalized_counts <- counts(dds, normalized = TRUE)
        
        # Get differential expression results
        res <- results(dds, contrast = c("condition", "tumor", "normal"))
        res_df <- as.data.frame(res)
        res_df$gene_id <- rownames(res_df)
        
        # Add gene names
        for (gene_id in res_df$gene_id) {
          if (gene_id %in% names(dataset$gene_names)) {
            gene_name <- dataset$gene_names[[gene_id]]
            
            # Extract statistics
            log2fc <- res_df[gene_id, "log2FoldChange"]
            p_value <- res_df[gene_id, "pvalue"]
            padj <- res_df[gene_id, "padj"]
            
            # Calculate means from normalized data
            normal_expr <- normalized_counts[gene_id, dataset$normal_samples]
            tumor_expr <- normalized_counts[gene_id, dataset$tumor_samples]
            
            normal_expr <- normal_expr[!is.na(normal_expr)]
            tumor_expr <- tumor_expr[!is.na(tumor_expr)]
            
            results <- bind_rows(results, tibble(
              project_id = project_id,
              gene_id = gene_id,
              gene_name = gene_name,
              mean_normal = mean(normal_expr),
              mean_tumor = mean(tumor_expr),
              median_normal = median(normal_expr),
              median_tumor = median(tumor_expr),
              log2fc = log2fc,
              p_value = p_value,
              padj = padj,
              n_normal = length(normal_expr),
              n_tumor = length(tumor_expr),
              comparison_type = "normal_vs_tumor"
            ))
          }
        }
        
        # Store normalized data for z-score calculation
        normalized_data[[project_id]] <- list(
          normalized_counts = normalized_counts,
          gene_names = dataset$gene_names,
          normal_samples = dataset$normal_samples,
          tumor_samples = dataset$tumor_samples
        )
        
      } else {
        # DESeq2 normalization without comparison (for z-score calculation)
        message("    Performing DESeq2 normalization for ", project_id, " (expression only)")
        
        # Create DESeq2 object with intercept-only design
        dds <- DESeqDataSetFromMatrix(
          countData = count_matrix,
          colData = sample_metadata,
          design = ~ 1
        )
        
        # Filter low count genes
        keep <- rowSums(counts(dds) >= 10) >= ceiling(0.1 * ncol(count_matrix))
        dds <- dds[keep, ]
        
        # Estimate size factors and dispersions
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersions(dds)
        
        # Get normalized counts
        normalized_counts <- counts(dds, normalized = TRUE)
        
        # Calculate expression statistics
        for (gene_id in rownames(normalized_counts)) {
          if (gene_id %in% names(dataset$gene_names)) {
            gene_name <- dataset$gene_names[[gene_id]]
            
            all_expr <- as.numeric(normalized_counts[gene_id, ])
            all_expr <- all_expr[!is.na(all_expr)]
            
            if (length(all_expr) >= 5) {
              results <- bind_rows(results, tibble(
                project_id = project_id,
                gene_id = gene_id,
                gene_name = gene_name,
                mean_expression = mean(all_expr),
                median_expression = median(all_expr),
                sd_expression = sd(all_expr),
                n_samples = length(all_expr),
                comparison_type = "expression_only"
              ))
            }
          }
        }
        
        # Store normalized data for z-score calculation
        normalized_data[[project_id]] <- list(
          normalized_counts = normalized_counts,
          gene_names = dataset$gene_names,
          normal_samples = dataset$normal_samples,
          tumor_samples = dataset$tumor_samples
        )
      }
      
    }, error = function(e) {
      message("    Error in DESeq2 analysis for ", project_id, ": ", e$message)
    })
  }
  
  return(list(results = results, normalized_data = normalized_data))
}

# Update the main analysis call
message("\n=== Performing DESeq2 normalization and differential expression analysis ===")
analysis_output <- perform_de_analysis(tcga_results)
analysis_results <- analysis_output$results
normalized_data <- analysis_output$normalized_data

# Save results
write_csv(analysis_results, file.path(output_dir, "tcga_chaperone_analysis_results.csv"))

# Create heatmaps
message("\\n=== Creating separate heatmaps ===")

# Helper function to get project names
get_project_name <- function(project_id) {
  project_name <- projects_info$name[projects_info$project_id == project_id]
  if (length(project_name) > 0) {
    # Shorten long names
    short_name <- str_replace(project_name, "^(.*?)\\s*-\\s*TCGA$", "\\1")
    return(short_name)
  } else {
    return(str_replace(project_id, "TCGA-", ""))
  }
}

# Create gene family annotation (common for both heatmaps)
gene_family_mapping <- setNames(gene_list$Family, gene_list$Name)

# 1. Z-score heatmap for all projects using DESeq2 normalized data
message("Creating z-score heatmap for all projects using DESeq2 normalized data...")

# Create tumor expression data from normalized counts
tumor_expression_data <- tibble()

for (project_id in names(normalized_data)) {
  project_data <- normalized_data[[project_id]]
  
  # Use tumor samples if available, otherwise all samples
  tumor_samples <- project_data$tumor_samples
  if (length(tumor_samples) == 0) {
    tumor_samples <- colnames(project_data$normalized_counts)
  }
  
  # Filter for available samples
  available_samples <- intersect(tumor_samples, colnames(project_data$normalized_counts))
  
  if (length(available_samples) >= 3) {
    for (gene_id in rownames(project_data$normalized_counts)) {
      if (gene_id %in% names(project_data$gene_names)) {
        gene_name <- project_data$gene_names[[gene_id]]
        
        tumor_expr <- as.numeric(project_data$normalized_counts[gene_id, available_samples])
        tumor_expr <- tumor_expr[!is.na(tumor_expr)]
        
        if (length(tumor_expr) >= 3) {
          tumor_expression_data <- bind_rows(tumor_expression_data, tibble(
            project_id = project_id,
            gene_name = gene_name,
            mean_tumor_expression = mean(tumor_expr)
          ))
        }
      }
    }
  }
}

# Create z-scores using DESeq2 normalized data
zscore_data <- tumor_expression_data %>%
  group_by(gene_name) %>%
  mutate(
    z_score = scale(log2(mean_tumor_expression + 1))[,1]
  ) %>%
  ungroup() %>%
  mutate(project_name = sapply(project_id, get_project_name))

# Create z-score matrix
zscore_matrix <- zscore_data %>%
  select(project_name, gene_name, z_score) %>%
  pivot_wider(names_from = gene_name, values_from = z_score) %>%
  column_to_rownames("project_name") %>%
  as.matrix()

# Create gene family annotation for z-score heatmap
matched_genes_zscore <- colnames(zscore_matrix)
gene_families_zscore <- sapply(matched_genes_zscore, function(gene) {
  if (gene %in% names(gene_family_mapping)) {
    return(gene_family_mapping[gene])
  } else {
    return("Unknown")
  }
})

col_anno_zscore <- data.frame(
  Family = gene_families_zscore,
  stringsAsFactors = FALSE
)
rownames(col_anno_zscore) <- matched_genes_zscore

# Create colors for families - sort families to ensure consistent color assignment across scripts
unique_families_zscore <- sort(unique(col_anno_zscore$Family))
# Use consistent family color assignment logic across all scripts
if (length(unique_families_zscore) <= 9) {
  family_colors_zscore <- setNames(brewer.pal(max(3, length(unique_families_zscore)), "Set1")[1:length(unique_families_zscore)], 
                                 unique_families_zscore)
} else {
  family_colors_zscore <- setNames(
    colorRampPalette(brewer.pal(9, "Set1"))(length(unique_families_zscore)),
    unique_families_zscore
  )
}

anno_colors_zscore <- list(Family = family_colors_zscore)

# Create z-score heatmap with same color scheme as log2FC
png(file.path(output_dir, "tcga_chaperone_zscore_heatmap.png"), 
    width = 16, height = 12, units = "in", res = 300)
pheatmap(
  zscore_matrix,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  annotation_col = col_anno_zscore,
  annotation_colors = anno_colors_zscore,
  main = "Chaperone Gene Expression Z-scores Across TCGA Projects",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = "45",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  border_color = NA,
  breaks = seq(-2, 2, length.out = 101),
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2.0", "-1.0", "0.0", "+1.0", "+2.0"),
  cellwidth = 25,
  cellheight = 15
)
dev.off()

# PDF version
pdf(file.path(output_dir, "tcga_chaperone_zscore_heatmap.pdf"), 
    width = 20, height = 14, onefile = FALSE)
print(pheatmap(
  zscore_matrix,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  annotation_col = col_anno_zscore,
  annotation_colors = anno_colors_zscore,
  main = "Chaperone Gene Expression Z-scores Across TCGA Projects",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = "45",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  border_color = NA,
  breaks = seq(-2, 2, length.out = 101),
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2.0", "-1.0", "0.0", "+1.0", "+2.0"),
  cellwidth = 25,
  cellheight = 15,
  silent = TRUE
))
dev.off()

message("Created z-score heatmap with ", nrow(zscore_matrix), " projects and ", ncol(zscore_matrix), " genes")
# Save CSV for z-score heatmap
message("Saving z-score heatmap CSV...")
zscore_df <- as.data.frame(zscore_matrix)
zscore_df$project_name <- rownames(zscore_df)
zscore_df <- zscore_df %>% select(project_name, everything())
write_csv(zscore_df, file.path(output_dir, "zscore_heatmap_data.csv"))
message("Saved: zscore_heatmap_data.csv")

# 2. Log2FC heatmap for projects with normal vs tumor comparison
message("Creating log2FC heatmap for projects with normal comparison...")

# Get normal vs tumor data (using the new analysis_results structure)
log2fc_data <- analysis_results %>%
  filter(comparison_type == "normal_vs_tumor") %>%
  select(project_id, gene_name, log2fc, p_value, padj) %>%
  filter(!is.na(log2fc), !is.infinite(log2fc)) %>%
  mutate(project_name = sapply(project_id, get_project_name))

if (nrow(log2fc_data) > 0) {
  # Create log2FC matrix
  log2fc_matrix <- log2fc_data %>%
    select(project_name, gene_name, log2fc) %>%
    pivot_wider(names_from = gene_name, values_from = log2fc) %>%
    column_to_rownames("project_name") %>%
    as.matrix()
  
  # Create significance annotation
  sig_data <- log2fc_data %>%
    select(project_name, gene_name, p_value) %>%
    pivot_wider(names_from = gene_name, values_from = p_value) %>%
    column_to_rownames("project_name") %>%
    as.matrix()
  
  # Create significance annotation matrix
  sig_annotation <- matrix("", nrow = nrow(log2fc_matrix), ncol = ncol(log2fc_matrix))
  rownames(sig_annotation) <- rownames(log2fc_matrix)
  colnames(sig_annotation) <- colnames(log2fc_matrix)
  
  # Fill in significance stars
  for (i in seq_len(nrow(log2fc_matrix))) {
    for (j in seq_len(ncol(log2fc_matrix))) {
      project_name <- rownames(log2fc_matrix)[i]
      gene_name <- colnames(log2fc_matrix)[j]
      if (project_name %in% rownames(sig_data) && gene_name %in% colnames(sig_data)) {
        p_val <- sig_data[project_name, gene_name]
        if (!is.na(p_val)) {
          if (p_val < 0.001) {
            sig_annotation[i, j] <- "***"
          } else if (p_val < 0.01) {
            sig_annotation[i, j] <- "**"
          } else if (p_val < 0.05) {
            sig_annotation[i, j] <- "*"
          }
        }
      }
    }
  }
  
  # Create gene family annotation for log2FC heatmap
  matched_genes_log2fc <- colnames(log2fc_matrix)
  gene_families_log2fc <- sapply(matched_genes_log2fc, function(gene) {
    if (gene %in% names(gene_family_mapping)) {
      return(gene_family_mapping[gene])
    } else {
      return("Unknown")
    }
  })
  
  col_anno_log2fc <- data.frame(
    Family = gene_families_log2fc,
    stringsAsFactors = FALSE
  )
  rownames(col_anno_log2fc) <- matched_genes_log2fc
  
  # Create colors for families - sort families to ensure consistent color assignment across scripts
  unique_families_log2fc <- sort(unique(col_anno_log2fc$Family))
  # Use consistent family color assignment logic across all scripts
  if (length(unique_families_log2fc) <= 9) {
    family_colors_log2fc <- setNames(brewer.pal(max(3, length(unique_families_log2fc)), "Set1")[1:length(unique_families_log2fc)], 
                                   unique_families_log2fc)
  } else {
    family_colors_log2fc <- setNames(
      colorRampPalette(brewer.pal(9, "Set1"))(length(unique_families_log2fc)),
      unique_families_log2fc
    )
  }
  
  anno_colors_log2fc <- list(Family = family_colors_log2fc)
  
  # Create log2FC heatmap
  png(file.path(output_dir, "tcga_chaperone_log2fc_heatmap.png"), 
      width = 16, height = 12, units = "in", res = 300)
  pheatmap(
    log2fc_matrix,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    display_numbers = sig_annotation,
    fontsize_number = 8,
    annotation_col = col_anno_log2fc,
    annotation_colors = anno_colors_log2fc,
    main = "Chaperone Gene Expression Changes (Tumor vs Normal)",
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 10,
    angle_col = "45",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    border_color = NA,
    breaks = seq(-2, 2, length.out = 101),
    legend_breaks = c(-2, -1, 0, 1, 2),
    legend_labels = c("-2.0", "-1.0", "0.0", "+1.0", "+2.0"),
    cellwidth = 25,
    cellheight = 15
  )
  dev.off()
  
  # PDF version
  pdf(file.path(output_dir, "tcga_chaperone_log2fc_heatmap.pdf"), 
      width = 20, height = 14, onefile = FALSE)
  print(pheatmap(
    log2fc_matrix,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    display_numbers = sig_annotation,
    fontsize_number = 8,
    annotation_col = col_anno_log2fc,
    annotation_colors = anno_colors_log2fc,
    main = "Chaperone Gene Expression Changes (Tumor vs Normal)",
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 10,
    angle_col = "45",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    border_color = NA,
    breaks = seq(-2, 2, length.out = 101),
    legend_breaks = c(-2, -1, 0, 1, 2),
    legend_labels = c("-2.0", "-1.0", "0.0", "+1.0", "+2.0"),
    cellwidth = 25,
    cellheight = 15,
    silent = TRUE
  ))
  dev.off()
  
  message("Created log2FC heatmap with ", nrow(log2fc_matrix), " projects and ", ncol(log2fc_matrix), " genes")
    # Save CSV for log2FC heatmap
  message("Saving log2FC heatmap CSV...")
  log2fc_df <- as.data.frame(log2fc_matrix)
  log2fc_df$project_name <- rownames(log2fc_df)
  log2fc_df <- log2fc_df %>% select(project_name, everything())
  write_csv(log2fc_df, file.path(output_dir, "log2fc_heatmap_data.csv"))
  message("Saved: log2fc_heatmap_data.csv")
} else {
  message("No projects with normal vs tumor comparison found")
}

# Create summary statistics
message("\\n=== Creating summary statistics ===")

# Summary for projects with normal vs tumor comparison
normal_tumor_data <- analysis_results %>%
  filter(comparison_type == "normal_vs_tumor")

if (nrow(normal_tumor_data) > 0) {
  normal_tumor_summary <- normal_tumor_data %>%
    group_by(gene_name) %>%
    summarise(
      projects_analyzed = n(),
      mean_log2fc = mean(log2fc, na.rm = TRUE),
      median_log2fc = median(log2fc, na.rm = TRUE),
      up_in_tumor = sum(log2fc > 0 & p_value < 0.05, na.rm = TRUE),
      down_in_tumor = sum(log2fc < 0 & p_value < 0.05, na.rm = TRUE),
      pct_up = round(100 * up_in_tumor / projects_analyzed, 1),
      pct_down = round(100 * down_in_tumor / projects_analyzed, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(pct_up))
  
  # Add gene family information
  normal_tumor_summary <- normal_tumor_summary %>%
    left_join(gene_list %>% select(Name, Family), by = c("gene_name" = "Name"))
  
  write_csv(normal_tumor_summary, file.path(output_dir, "normal_tumor_summary.csv"))
}

# Summary for all projects
overall_summary <- tibble(
  project_id = names(normalized_data),
  total_samples = map_int(normalized_data, ~ ncol(.x$normalized_counts)),
  normal_samples = map_int(normalized_data, ~ length(.x$normal_samples)),
  tumor_samples = map_int(normalized_data, ~ length(.x$tumor_samples)),
  chaperone_genes_found = map_int(normalized_data, ~ nrow(.x$normalized_counts)),
  has_normal_comparison = map_lgl(normalized_data, ~ length(.x$normal_samples) >= 3 && length(.x$tumor_samples) >= 3)
)

write_csv(overall_summary, file.path(output_dir, "project_summary.csv"))

# Print summary
message("\\n=== Analysis Summary ===")
message("Total TCGA projects processed: ", length(normalized_data))
message("Projects with normal vs tumor comparison: ", sum(overall_summary$has_normal_comparison))
message("Projects with expression data only: ", sum(!overall_summary$has_normal_comparison))
message("Average chaperone genes found per project: ", round(mean(overall_summary$chaperone_genes_found), 1))

if (exists("normal_tumor_summary")) {
  message("\\nTop upregulated chaperones in tumors:")
  print(head(normal_tumor_summary %>% arrange(desc(pct_up)) %>% select(gene_name, pct_up, mean_log2fc), 5))
  
  message("\\nTop downregulated chaperones in tumors:")
  print(head(normal_tumor_summary %>% arrange(desc(pct_down)) %>% select(gene_name, pct_down, mean_log2fc), 5))
}

message("\\nAnalysis complete! Results saved to: ", output_dir)
message("Completed at: ", Sys.time())

# Close logging
sink()
close(con)
