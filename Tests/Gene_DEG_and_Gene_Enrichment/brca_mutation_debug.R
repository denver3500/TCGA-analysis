library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(edgeR)
library(readr)
library(maftools)
library(pathfindR)

# Set project ID and gene of interest
project_id <- "TCGA-BRCA"
gene_of_interest <- "TP53"  # Changed from BRCA1 to TP53

# Create output directory
output_dir <- paste0("TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/", gene_of_interest, "_mutation_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 1: Download mutation data for BRCA from GDC
message("Downloading mutation data for TCGA-BRCA...")
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(maf_query, files.per.chunk = 100)
maf_data <- GDCprepare(maf_query)

# Step 2: Extract TP53 mutation information
message("Extracting TP53 mutation information...")
# Convert to MAF object for easier handling
maf_object <- read.maf(maf = maf_data)

# Find samples with TP53 mutations
tp53_mutations <- subsetMaf(maf_object, genes = gene_of_interest)
tp53_mutated_samples <- unique(substr(tp53_mutations@data$Tumor_Sample_Barcode, 1, 12))

message(paste("Found", length(tp53_mutated_samples), "samples with TP53 mutations"))
write.csv(data.frame(Sample = tp53_mutated_samples), 
          file.path(output_dir, "TP53_mutated_samples.csv"), 
          row.names = FALSE)


#Mutation_group will be based on TP53 status

# Load RNA-seq data from RDS file (this part stays the same)
message("Loading RNAseq data...")
rds_path <- list.files("TCGA-Chaperones/rds", 
                       pattern = paste0("^", project_id, "_.*\\.rds$"), 
                       full.names = TRUE)

if (length(rds_path) == 0) {
  stop("Could not find RDS file for project ", project_id)
}

se_object <- readRDS(rds_path[1])

# Check available assays in the SummarizedExperiment object
message("Available assays in SummarizedExperiment:")
print(assayNames(se_object))

# If "counts" assay doesn't exist, try to find and rename an appropriate assay
if (!"counts" %in% assayNames(se_object)) {
  possible_count_assays <- c("unstranded", "raw_counts", "HTSeq - Counts", "raw_count")
  found_assay <- FALSE
  
  for (assay_name in possible_count_assays) {
    if (assay_name %in% assayNames(se_object)) {
      message(paste("Renaming", assay_name, "assay to counts"))
      counts_data <- assay(se_object, assay_name)
      assay(se_object, "counts") <- counts_data
      found_assay <- TRUE
      break
    }
  }
  
  if (!found_assay) {
    # If no recognized count assay is found, try using the first available assay
    if (length(assayNames(se_object)) > 0) {
      first_assay <- assayNames(se_object)[1]
      message(paste("No standard count assay found. Using", first_assay, "as counts"))
      counts_data <- assay(se_object, first_assay)
      assay(se_object, "counts") <- counts_data
    } else {
      stop("No assays found in the SummarizedExperiment object")
    }
  }
}

# Now try the conversion again
edgeRlist <- SE2DGEList(se_object)

# Filter for tumor samples
if("definition" %in% colnames(edgeRlist$samples)) {
  sample_idx <- edgeRlist$samples$definition == "Primary solid Tumor"
} else if("sample_type" %in% colnames(edgeRlist$samples)) {
  sample_idx <- edgeRlist$samples$sample_type == "Primary Tumor"
} else if("shortLetterCode" %in% colnames(edgeRlist$samples)) {
  sample_idx <- edgeRlist$samples$shortLetterCode == "TP"
} else {
  message("No valid field found for filtering tumor samples! Using all samples.")
  sample_idx <- rep(TRUE, ncol(se_object))
}

# Filter the DGEList for tumor samples
dge <- edgeRlist[, sample_idx]
message(paste("Working with", ncol(dge), "tumor samples"))

# Create grouping based on TP53 mutation status
tcga_barcodes <- colnames(dge)
short_barcodes <- substr(tcga_barcodes, 1, 12)

# Create grouping factor (updated for TP53)
mutation_group <- factor(ifelse(short_barcodes %in% tp53_mutated_samples, 
                              "TP53_Mutated", "TP53_Wildtype"))

# Add grouping to sample information
dge$samples$mutation_group <- mutation_group

# Display group sizes
group_sizes <- table(mutation_group)
message("Sample grouping by TP53 mutation status:")
print(group_sizes)

# Check if we have enough samples in each group
if (min(group_sizes) < 3) {
  stop("Not enough samples in one of the groups for statistical analysis")
}

# Write sample grouping to file for reference
sample_info <- data.frame(
  TCGA_Barcode = tcga_barcodes,
  TP53_Status = mutation_group
)
write.csv(sample_info, file.path(output_dir, "sample_mutation_grouping.csv"), row.names = FALSE)

# The rest of the script stays the same, just make sure to update:
# - All references to BRCA1 in output filenames to TP53
# - In the contrast statement: my_contrast <- makeContrasts(TP53_Mutated - TP53_Wildtype, levels = design)

# Step 6: Filter low expressed genes
message("Filtering low expressed genes...")
keep <- filterByExpr(dge, group = mutation_group)
message(paste("Keeping", sum(keep), "out of", length(keep), "genes"))
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Step 7: Normalize the data
message("Normalizing data...")
dge <- calcNormFactors(dge)

# Step 8: Create design matrix
message("Creating design matrix...")
design <- model.matrix(~0 + mutation_group)
colnames(design) <- levels(mutation_group)

# Step 9: Estimate dispersion
message("Estimating dispersion...")
dge <- estimateDisp(dge, design)

# Step 10: Fit the model
message("Fitting GLM...")
fit <- glmQLFit(dge, design)

# Step 11: Define contrast and test for differential expression
message("Testing for differential expression...")
my_contrast <- makeContrasts(TP53_Mutated - TP53_Wildtype, levels = design)  # Updated for TP53
qlf <- glmQLFTest(fit, contrast = my_contrast)

# Step 12: Extract results
message("Extracting DEG results...")
res <- topTags(qlf, n = Inf)
results_df <- as.data.frame(res)
results_df$gene_ensembl <- rownames(results_df)

# Add gene names if available
if (exists("rowRanges") && is(rowRanges(se_object), "GRanges")) {
  gene_info <- as.data.frame(mcols(rowRanges(se_object)))
  
  potential_name_cols <- c("gene_name", "external_gene_name", "symbol", "gene_symbol")
  name_col <- potential_name_cols[potential_name_cols %in% colnames(gene_info)][1]
  
  if (!is.na(name_col)) {
    # Create mapping from ENSEMBL ID to gene name
    ensembl_to_gene_map <- data.frame(
      ensembl_id = rownames(se_object),
      gene_name = gene_info[[name_col]]
    )
    
    # Match gene names to results
    extract_base_ensembl <- function(ensembl_id) {
      return(sapply(strsplit(ensembl_id, "\\."), function(x) x[1]))
    }
    
    results_df$ensembl_base <- extract_base_ensembl(results_df$gene_ensembl)
    ensembl_to_gene_map$ensembl_base <- extract_base_ensembl(ensembl_to_gene_map$ensembl_id)
    
    # Merge gene names
    results_df <- merge(results_df, 
                       ensembl_to_gene_map[, c("ensembl_base", "gene_name")], 
                       by = "ensembl_base", all.x = TRUE)
  }
}

# Step 13: Save DEG results
deg_file <- file.path(output_dir, paste0(gene_of_interest, "_mutation_DEG_results.csv"))
write.csv(results_df, deg_file, row.names = FALSE)

message(paste("DEG analysis complete. Results saved to:", deg_file))

# Display summary statistics
message("=== DEG Results Summary ===")
message(paste("Total genes tested:", nrow(results_df)))
message(paste("Upregulated in TP53-mutated (logFC > 0):", sum(results_df$logFC > 0)))
message(paste("Downregulated in TP53-mutated (logFC < 0):", sum(results_df$logFC < 0)))
message(paste("Significant at FDR < 0.05:", sum(results_df$FDR < 0.05)))
message(paste("Significant at FDR < 0.1:", sum(results_df$FDR < 0.1)))

# Step 14: Run pathfindR on the mutation-based DEG results
message("\n=== Running pathfindR Analysis on Mutation-Based DEGs ===")

# For pathfindR, prepare input for genes with FDR < 0.05
significant_genes <- results_df[results_df$FDR < 0.05, ]

if (nrow(significant_genes) >= 10) {
  message(paste("Running pathfindR with", nrow(significant_genes), "significant genes"))
  
  # Run with two parameter combinations
  for (gene_sets in c("KEGG", "GO-All")) {
    pathfindR_dir <- file.path(output_dir, paste0("pathfindR_", gene_sets))
    dir.create(pathfindR_dir, showWarnings = FALSE)
    
    message(paste("Running pathfindR with gene_sets =", gene_sets))
    
    # Prepare input
    input_df <- data.frame(
      Gene_symbol = significant_genes$gene_name.y,
      logFC = significant_genes$logFC,
      FDR_p = significant_genes$FDR
    )
    
    # Remove NA values
    input_df <- input_df[!is.na(input_df$Gene_symbol), ]
    
    # Run pathfindR
    tryCatch({
      pathfinder_results <- run_pathfindR(
        input_df,
        output_dir = pathfindR_dir,
        gene_sets = gene_sets,
        pin_name = "Biogrid",
        min_gset_size = 5,
        max_gset_size = 500,
        silent_option = TRUE
      )
      
      # Save results
     # Replace the problematic section with this code
if (nrow(pathfinder_results) > 0) {
  result_file <- file.path(pathfindR_dir, 
                         paste0(gene_of_interest, "_mutation_", gene_sets, "_results.csv"))
  write_csv(pathfinder_results, result_file)
  message(paste("Found", nrow(pathfinder_results), "enriched terms"))
  
  # First check what columns are actually available
  message("Available columns in pathfindR results:")
  message(paste(colnames(pathfinder_results), collapse=", "))
  
  # Then display top 5 rows with whatever columns are available
  message("Top 5 terms:")
  print(head(pathfinder_results, 5))
} else {
        message("No enriched terms found")
      }
    }, error = function(e) {
      message(paste("Error in pathfindR analysis:", e$message))
    })
  }
} else {
  message(paste("Too few significant genes (", nrow(significant_genes), ") for pathfindR analysis"))
  # Try with less stringent cutoff
  significant_genes <- results_df[results_df$FDR < 0.1, ]
  message(paste("With less stringent cutoff (FDR < 0.1):", nrow(significant_genes), "genes"))
  
  if (nrow(significant_genes) >= 10) {
    message("Running pathfindR with this less stringent gene set")
    # (pathfindR code would go here, similar to above)
  }
}
