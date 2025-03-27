library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(corrplot)

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
  # Extract the project ID from filenames like "TCGA-BRCA_transcriptomic_exp.rds"
  # This splits at the underscore and takes the first part
  filename <- basename(file)
  project_id <- strsplit(filename, "_")[[1]][1]
  
  return(project_id %in% filtered_projects$project_id)
})]

# Print summary of filtered files
message(paste("Found", length(filtered_rds_files), "matching RDS files out of", length(rds_files), "total"))

for(file in filtered_rds_files) {
  # Extract just the project ID from the filename for use in output files
  filename <- basename(file)
  project <- strsplit(filename, "_")[[1]][1]
  
  cat("Processing file:", file, "\n")
  # Read in the expression data - use the full path directly
  brca.exp <- readRDS(file)
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
  
  # Get gene names from the gene_list CSV instead of hardcoding
  gene_names <- gene_list$Name
  
  # Get gene IDs for these genes
  gene_ids <- gene_metadata[gene_metadata[, "gene_name"] %in% gene_names, "gene_id"]
  
  # Check if any gene IDs were found
  if(length(gene_ids) == 0) {
    message(paste("Warning: No matching gene IDs found for project", project))
    next
  }
  # Preprocess data
 
  dataNormFactors <- calcNormFactors(edgeRlist.subset)
  dataCPM <- cpm(dataNormFactors)
  
  # Subset CPM data for selected genes
  dataCPM_sub <- dataCPM[rownames(dataCPM) %in% gene_ids, ]
  
 # Compute the pairwise correlation matrix
  cor_matrix <- cor(t(dataCPM_sub), method = "pearson")
  
  # Use the gene_list.csv which you already have loaded instead of another mapping file
  # This assumes your gene_list.csv has the information you need
  gene_id_to_name <- setNames(gene_metadata[, "gene_name"], gene_metadata[, "gene_id"])
  
  # Map the rownames (gene_ids) to gene names
  gene_names_for_plot <- sapply(rownames(dataCPM_sub), function(id) {
    if (id %in% names(gene_id_to_name)) {
      return(gene_id_to_name[id])
    } else {
      return(id) # If no mapping found, keep the original
    }
  })
  
  # Set row and column names
  rownames(cor_matrix) <- gene_names_for_plot
  colnames(cor_matrix) <- gene_names_for_plot
 # Generate the correlation plot and save as PNG
  file_name <- file.path(
    "TCGA-Expression",
    "pictures",
    "correlations",
    paste0("corrplot_", project, ".png"))
  
  # Create directory if it doesn't exist
  dir.create(file.path("TCGA-Expression", "pictures", "correlations"), showWarnings = FALSE, recursive = TRUE)
  
  png(filename = file_name, width = 1000, height = 800)
  corrplot(cor_matrix, method = "color", type = "lower",  title = paste("Correlation Comparison for", project),  mar = c(0, 0, 1, 0), 
           tl.col = "black", tl.srt = 45)
  dev.off()
  
 # Generate the correlation plot and save as PNG
  file_name <- file.path(
    "TCGA-Chaperones",
    "pictures",
    "correlations",
    paste0("corrplot_", project, ".png"))
  
  # Create directory if it doesn't exist
  dir.create(file.path("TCGA-Chaperones", "pictures", "correlations"), showWarnings = FALSE, recursive = TRUE)
  
  png(filename = file_name, width = 1000, height = 800)
  corrplot(cor_matrix, method = "color", type = "lower",  title = paste("Correlation Comparison for", project),  mar = c(0, 0, 1, 0), 
           tl.col = "black", tl.srt = 45)
  dev.off()
  
  # Generate the correlation plot and save as PNG
  file_name <- file.path(
    "TCGA-Chaperones",
    "pictures",
    "correlations_aoe",
    paste0("corrplot_", project, "_AOE", ".png"))
  
  # Create directory if it doesn't exist
  dir.create(file.path("TCGA-Chaperones", "pictures", "correlations_aoe"), showWarnings = FALSE, recursive = TRUE)
  
  png(filename = file_name, width = 1000, height = 800)
  corrplot(cor_matrix, method = "color", type = "lower", order = "AOE", title = paste("Correlation Comparison for", project),  mar = c(0, 0, 1, 0),
           tl.col = "black", tl.srt = 45)
  dev.off()
}