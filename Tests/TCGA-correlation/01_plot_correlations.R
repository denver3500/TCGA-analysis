library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(corrplot)

# List all .rds files in the working directory
rds_files <- list.files("raw_data/", pattern = "*.rds") 

for(file in rds_files) {
  # Remove the .rds extension to use as the project name
  project <- sub("\\.rds$", "", file)
  cat("Processing file:", file, "\n")
  # Read in the expression data
  brca.exp <- readRDS(file.path("raw_data", file))
  assayNames(brca.exp)[assayNames(brca.exp) == "unstranded"] <- "counts"
  edgeRlist <- SE2DGEList(brca.exp)
  
  # Check if each field exists in edgeRlist$samples and create a condition accordingly. Omg, trying to filter with multiple conditions took me 3 hours. People should use same annotation for each project >_<
  if("definition" %in% colnames(edgeRlist$samples)) {
    message("Using 'definition' field")
    sample.idx <- edgeRlist$samples$definition == "Primary solid Tumor"
  } else if("sample_type" %in% colnames(edgeRlist$samples)) {
    message("Using 'sample_type' field")
    sample.idx <- edgeRlist$samples$sample_type == "Primary Tumor"
  } else if("shortLetterCode" %in% colnames(edgeRlist$samples)) {
    message("Using 'shortLetterCode' field")
    sample.idx <- edgeRlist$samples$shortLetterCode == "TP"
  } else {
    message("No valid field found in edgeRlist$samples for filtering!")
    next
  }
  
  # Use the vector to subset the DGEList
  edgeRlist.subset <- edgeRlist[, sample.idx]
  
  # Extract gene metadata and subset data for TP samples
  gene_metadata <- as.matrix(mcols(rowRanges(brca.exp))[, c("gene_id", "gene_name")])
  
  # Gene selection
  gene_names <- c("HSP90AA1", "HSPA1A", "STIP1", "NR3C1", "PTGES3", 
                  "FKBP5", "FKBP4", "HSP90AB1", "CDC37", "CDK4")
  gene_ids <- gene_metadata[gene_metadata[, "gene_name"] %in% gene_names, "gene_id"]
  
  # Preprocess data
 
  dataNormFactors <- calcNormFactors(edgeRlist.subset)
  dataCPM <- cpm(dataNormFactors)
  
  # Subset CPM data for selected genes
  dataCPM_sub <- dataCPM[rownames(dataCPM) %in% gene_ids, ]
  
  # Compute the pairwise correlation matrix
  cor_matrix <- cor(t(dataCPM_sub), method = "pearson")
  mapping <- read.csv("TCGA-correlation/gene_protein_name.csv", stringsAsFactors = FALSE)
  # Remove extra spaces and standardize names to lowercase for matching
  mapping[, 1] <- tolower(trimws(mapping[, 1]))
  gene_names <- tolower(trimws(gene_names))
  gene_names <- mapping[match(gene_names, mapping[, 1]), 2]
  rownames(cor_matrix) <- gene_names
  colnames(cor_matrix) <- gene_names
  # Generate the correlation plot and save as PNG
  file_name <- file.path(
    "pictures",
    "correlations",
    paste0("corrplot_", project, ".png"))
  png(filename = file_name, width = 1000, height = 800)
  corrplot(cor_matrix, method = "color", type = "lower",  title = paste("Correlation Comparison for", project),  mar = c(0, 0, 1, 0), 
           tl.col = "black", tl.srt = 45)
  dev.off()
  
  # Generate the correlation plot and save as PNG
  file_name <- file.path(
    "pictures",
    "correlations_AOE",
    paste0("corrplot_", project, "_AOE", ".png"))
  png(filename = file_name, width = 1000, height = 800)
  corrplot(cor_matrix, method = "color", type = "lower", order = "AOE", title = paste("Correlation Comparison for", project),  mar = c(0, 0, 1, 0),
           tl.col = "black", tl.srt = 45)
  dev.off()
}