library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(dplyr)

# Set up logging to both console and file
log_file <- "TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/debug_HSP90AA1_BRCA.txt"
con <- file(log_file, "w")
sink(con, split=TRUE)  # Write output to both console and file

# Define the specific project and gene to analyze
project_id <- "TCGA-BRCA"
gene_name <- "BRCA1"
gene_ensembl <- "ENSG00000012048"

message("====== DEBUG ANALYSIS START ======")
message(paste("Project:", project_id))
message(paste("Gene:", gene_name, "(", gene_ensembl, ")"))

# Find and load the RDS file
rds_path <- list.files("TCGA-Chaperones/rds", 
                      pattern = paste0("^", project_id, "_.*\\.rds$"), 
                      full.names = TRUE)

if (length(rds_path) == 0) {
  stop(paste("Could not find RDS file for project", project_id))
}

message(paste("Loading RDS file:", basename(rds_path[1])))
se_object <- readRDS(rds_path[1])

# Print general information about the dataset
message("=== Dataset Information ===")
message(paste("Dimensions:", nrow(se_object), "genes,", ncol(se_object), "samples"))
message(paste("Assay names:", paste(assayNames(se_object), collapse=", ")))

# Check and rename assays if necessary
if ("unstranded" %in% assayNames(se_object)) {
  message("Renaming 'unstranded' assay to 'counts'")
  assayNames(se_object)[assayNames(se_object) == "unstranded"] <- "counts"
}

# Extract metadata columns
message("=== Sample Metadata Columns ===")
edgeRlist <- tryCatch({
  SE2DGEList(se_object)
}, error = function(e) {
  message(paste("Error converting to DGEList:", e$message))
  return(NULL)
})

if (!is.null(edgeRlist)) {
  message(paste("Metadata columns:", paste(colnames(edgeRlist$samples), collapse=", ")))
  
  # Check for tumor/normal annotation columns
  potential_cols <- c("definition", "sample_type", "shortLetterCode")
  for (col in potential_cols) {
    if (col %in% colnames(edgeRlist$samples)) {
      message(paste("Values in", col, ":"))
      print(table(edgeRlist$samples[[col]]))
    }
  }
  
  # Filter for tumor samples
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
    message("No valid field found for filtering tumor samples! Using all samples.")
    sample_idx <- rep(TRUE, ncol(se_object))
  }
  
  message(paste("Number of tumor samples:", sum(sample_idx), "out of", length(sample_idx), "total"))
  
  # Filter the SummarizedExperiment
  se_object <- se_object[, sample_idx]
  message(paste("Filtered to", ncol(se_object), "tumor samples"))
} else {
  message("Could not access sample metadata - using all samples")
}

# Check for gene in the dataset
message("=== Gene Lookup ===")
assay_data <- assay(se_object)

# First try exact match
gene_found <- gene_ensembl %in% rownames(assay_data)
message(paste("Exact gene match found:", gene_found))

if (!gene_found) {
  # Look for genes that start with the ENSEMBL ID followed by a period
  pattern <- paste0("^", gene_ensembl, "\\.[0-9]+$")
  matching_genes <- grep(pattern, rownames(assay_data), value = TRUE)
  
  if (length(matching_genes) > 0) {
    message(paste("Found gene with version:", matching_genes[1]))
    gene_ensembl <- matching_genes[1]  # Use the first matching gene
    gene_found <- TRUE
  } else {
    stop(paste("Gene", gene_ensembl, "(", gene_name, ") not found in", project_id))
  }
}

# Extract gene expression
message("=== Gene Expression Analysis ===")
gene_expr <- assay_data[gene_ensembl, ]

# Print expression summary
message("Expression summary statistics:")
print(summary(gene_expr))

# Calculate Z-score and group samples
gene_zscore <- scale(gene_expr)
message("Z-score summary statistics:")
print(summary(gene_zscore))

# Create grouping based on Z-score
group_factor <- factor(ifelse(gene_zscore > 0, "High", "Low"))
group_table <- table(group_factor)
message("Sample grouping by Z-score (count):")
print(group_table)
message("Sample grouping by Z-score (percentage):")
print(prop.table(group_table) * 100)

# Print detailed expression by group
message("=== Expression Distribution by Group ===")
group_data <- data.frame(
  Expression = gene_expr,
  Z_Score = gene_zscore,
  Group = group_factor
)

message("High expression group summary:")
print(summary(group_data$Expression[group_data$Group == "High"]))
message("Low expression group summary:")
print(summary(group_data$Expression[group_data$Group == "Low"]))

# Export sample IDs with their groups
sample_groups <- data.frame(
  Sample_ID = colnames(se_object),
  Expression = gene_expr,
  Z_Score = gene_zscore,
  Group = group_factor
)
write.csv(sample_groups, paste0("TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/", project_id, "_", gene_name, "_sample_groups.csv"), row.names = FALSE)
message(paste("Exported sample grouping to CSV file"))

# Perform differential expression analysis
message("=== Differential Expression Analysis ===")
message("Creating DGEList object...")
dge <- DGEList(counts = assay_data)

# Filter low expressed genes
message("Filtering low expressed genes...")
keep <- filterByExpr(dge, group = group_factor)
message(paste("Keeping", sum(keep), "out of", length(keep), "genes (", round(sum(keep)/length(keep)*100, 1), "%)"))
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize the data
message("Normalizing data...")
dge <- calcNormFactors(dge)

# Create design matrix
message("Creating design matrix...")
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)
message("Design matrix dimensions:")
print(dim(design))
message("Design matrix column names:")
print(colnames(design))

# Estimate dispersion
message("Estimating dispersion...")
dge <- estimateDisp(dge, design)
message("BCV summary:")
print(summary(sqrt(dge$tagwise.dispersion)))

# Fit the model
message("Fitting GLM...")
fit <- glmQLFit(dge, design)

# Define contrast
message("Creating contrast matrix...")
my_contrast <- makeContrasts(High - Low, levels = design)
print(my_contrast)

# Test for differential expression
message("Testing for differential expression...")
qlf <- glmQLFTest(fit, contrast = my_contrast)

# Extract results
message("Extracting top differentially expressed genes...")
res <- topTags(qlf, n = 20)
print(res)

# Add gene info
results_df <- as.data.frame(topTags(qlf, n = Inf))
results_df$gene_ensembl <- rownames(results_df)
results_df$grouping_gene <- gene_name
results_df$project_id <- project_id

# Save full results
output_file <- paste0("TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/", project_id, "_", gene_name, "_DEG_results_debug.csv")
write.csv(results_df, output_file, row.names = FALSE)
message(paste("Saved full results to:", output_file))

# Summary statistics about results
message("=== DEG Results Summary ===")
message(paste("Total genes tested:", nrow(results_df)))
message(paste("Genes with FDR < 0.05:", sum(results_df$FDR < 0.05)))
message(paste("Genes with FDR < 0.1:", sum(results_df$FDR < 0.1)))
message(paste("Upregulated genes (logFC > 0):", sum(results_df$logFC > 0)))
message(paste("Downregulated genes (logFC < 0):", sum(results_df$logFC < 0)))

message("====== DEBUG ANALYSIS COMPLETE ======")

# Close the log connection
sink() 
close(con)

message("Debug log saved to:", log_file)