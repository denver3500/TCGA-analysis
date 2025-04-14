library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(dplyr)

# Create output directory for results
output_dir <- "TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/DEG_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the CSV file with projects and their gene counts
project_names <- read.csv("TCGA-Chaperones/significant_higher_in_tumor_gene_counts.csv", stringsAsFactors = FALSE)
gene_list <- read.csv("TCGA-Chaperones/gene_list.csv", stringsAsFactors = FALSE)

# Filter projects with 10 or more genes significantly higher in tumor
filtered_projects <- project_names[project_names$higher_in_tumor_count >= 10, ]

# Print summary of filtered projects
message(paste("Total projects with 10+ genes significantly higher in tumor:", nrow(filtered_projects)))
message("Projects included:")
print(filtered_projects$project_id)

# Filter rds files to only include those projects
rds_files <- list.files("TCGA-Chaperones/rds", pattern = "*.rds", full.names = TRUE)

filtered_rds_files <- rds_files[sapply(rds_files, function(file) {
  filename <- basename(file)
  project_id <- strsplit(filename, "_")[[1]][1]
  return(project_id %in% filtered_projects$project_id)
})]

# Function to perform DEG analysis between high and low gene expression groups
perform_deg_analysis <- function(se_object, gene_ensembl, gene_name, project_id) {
  # Extract gene expression data
  assay_data <- assay(se_object)
  
  # Check if gene exists in dataset (accounting for version numbers)
  # First try exact match
  gene_found <- gene_ensembl %in% rownames(assay_data)
  
  # If not found, try matching with version numbers
  if (!gene_found) {
    # Look for genes that start with the ENSEMBL ID followed by a period
    pattern <- paste0("^", gene_ensembl, "\\.[0-9]+$")
    matching_genes <- grep(pattern, rownames(assay_data), value = TRUE)
    
    if (length(matching_genes) > 0) {
      message(paste("Found gene", gene_ensembl, "with version:", matching_genes[1]))
      gene_ensembl <- matching_genes[1]  # Use the first matching gene
      gene_found <- TRUE
    } else {
      message(paste("Gene", gene_ensembl, "(", gene_name, ") not found in", project_id))
      return(NULL)
    }
  }
  
  # Get expression values for the target gene
  gene_expr <- assay_data[gene_ensembl, ]
  
  # Calculate Z-score for the gene across all samples
  gene_zscore <- scale(gene_expr)
  
  # Divide samples into high and low expression groups based on Z-score
  group_factor <- factor(ifelse(gene_zscore > 0, "High", "Low"))
  
  # Check if we have enough samples in each group
  table_groups <- table(group_factor)
  if (any(table_groups < 3)) {
    message(paste("Not enough samples in each group for gene", gene_name, "in", project_id))
    return(NULL)
  }
  
  # Create a DGEList object for edgeR
  dge <- DGEList(counts = assay_data)
  
  # Filter low expressed genes
  keep <- filterByExpr(dge, group = group_factor)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Normalize the data
  dge <- calcNormFactors(dge)
  
  # Create design matrix
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)
  
  # Estimate dispersion
  dge <- estimateDisp(dge, design)
  
  # Fit the model
  fit <- glmQLFit(dge, design)
  
  # Define contrast for High vs Low
  my_contrast <- makeContrasts(High - Low, levels = design)
  
  # Test for differential expression
  qlf <- glmQLFTest(fit, contrast = my_contrast)
  
  # Extract results
  res <- topTags(qlf, n = Inf)
  
  # Add gene info and return
  results_df <- as.data.frame(res)
  results_df$gene_ensembl <- rownames(results_df)
  results_df$grouping_gene <- gene_name
  results_df$project_id <- project_id
  
  return(results_df)
}

# Process each project
for (file_path in filtered_rds_files) {
  # Extract project ID from filename
  filename <- basename(file_path)
  project_id <- strsplit(filename, "_")[[1]][1]
  
  message(paste("Processing", project_id))
  
  # Load the RDS file
  se_object <- readRDS(file_path)
  
  # Check and rename assays if necessary for SE2DGEList compatibility
  if ("unstranded" %in% assayNames(se_object)) {
    message("Renaming 'unstranded' assay to 'counts'")
    assayNames(se_object)[assayNames(se_object) == "unstranded"] <- "counts"
  }
  
  # Filter for tumor samples only
  # Convert to edgeR object temporarily to access sample annotations
  tryCatch({
    edgeRlist <- SE2DGEList(se_object)
    
    # Check for different field names and filter accordingly
    if("definition" %in% colnames(edgeRlist$samples)) {
      message("Using 'definition' field for tumor filtering")
      sample_idx <- edgeRlist$samples$definition == "Primary solid Tumor"
    } else if("sample_type" %in% colnames(edgeRlist$samples)) {
      message("Using 'sample_type' field for tumor filtering")
      sample_idx <- edgeRlist$samples$sample_type == "Primary Tumor"
    } else if("shortLetterCode" %in% colnames(edgeRlist$samples)) {
      message("Using 'shortLetterCode' field for tumor filtering")
      sample_idx <- edgeRlist$samples$shortLetterCode == "TP"
    } else {
      message("No valid field found for filtering tumor samples in", project_id, "! Using all samples.")
      sample_idx <- rep(TRUE, ncol(se_object))
    }
    
    # Check if we have enough tumor samples
    if(sum(sample_idx) < 10) {
      message(paste("Too few tumor samples (", sum(sample_idx), ") in", project_id, "- skipping."))
      next
    }
    
  }, error = function(e) {
    message(paste("Error converting SE to DGEList:", e$message))
    message("Using all samples without filtering.")
    sample_idx <- rep(TRUE, ncol(se_object))
  })
  
  # Filter the SummarizedExperiment to include only tumor samples
  se_object <- se_object[, sample_idx]
  message(paste("  Filtered to", ncol(se_object), "tumor samples"))
  
  # Create a directory for this project
  project_dir <- file.path(output_dir, project_id)
  dir.create(project_dir, showWarnings = FALSE)
  
  # Process each gene in the gene list
  for (i in 1:nrow(gene_list)) {
    gene_ensembl <- gene_list$ENSEMBL[i]
    gene_name <- gene_list$Name[i]
    
    message(paste("  Analyzing gene:", gene_name))
    
    # Perform differential expression analysis
    tryCatch({
      deg_results <- perform_deg_analysis(se_object, gene_ensembl, gene_name, project_id)
      
      # Save results if not NULL
      if (!is.null(deg_results)) {
        output_file <- file.path(project_dir, paste0(gene_name, "_DEG_results.csv"))
        write.csv(deg_results, output_file, row.names = FALSE)
        message(paste("    Results saved to:", output_file))
      }
    }, error = function(e) {
      message(paste("Error processing gene", gene_name, "in", project_id, ":", e$message))
    })
  }
  
  message(paste("Completed analysis for", project_id))
}

message("Analysis complete!")