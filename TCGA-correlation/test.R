library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(corrplot)

# List all .rds files in the working directory
rds_files <- list.files(pattern = "\\.rds$")

for(file in rds_files) {
  # Remove the .rds extension to use as the project name
  project <- sub("\\.rds$", "", file)
  cat("Processing file:", file, "\n")
  
  # Read in the expression data
  brca.exp <- readRDS(file)
  
  # Check if there are any TP samples; if not, skip the file
  if(sum(brca.exp$shortLetterCode == "TP") == 0) {
    cat("No TP samples found in", file, "- skipping.\n")
    next
  }
  
  # Extract gene metadata and subset data for TP samples
  gene_metadata <- as.matrix(mcols(rowRanges(brca.exp))[, c("gene_id", "gene_name")])
  brca.exp.TP <- brca.exp[, brca.exp$shortLetterCode == "TP"]
  
  # Add a dummy 'Disease' column if necessary
  if (!"Disease" %in% colnames(colData(brca.exp.TP))) {
    colData(brca.exp.TP)$Disease <- rep("unknown", ncol(brca.exp.TP))
  }
  
  # Gene selection
  gene_names <- c("HSP90AA1", "HSPA1A", "STIP1", "NR3C1", "PTGES3", 
                  "FKBP5", "FKBP4", "HSP90AB1", "CDC37", "CDK4")
  gene_ids <- gene_metadata[gene_metadata[, "gene_name"] %in% gene_names, "gene_id"]
  
  # Preprocess data
  dataPrep <- TCGAanalyze_Preprocessing(object = brca.exp.TP, cor.cut = 0.6)
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataPrep, method = "quantile", qnt.cut = 0.25)
  dataNormFactors <- calcNormFactors(dataFilt)
  dataCPM <- cpm(dataFilt)
  
  # Subset CPM data for selected genes
  dataCPM_sub <- dataCPM[rownames(dataCPM) %in% gene_ids, ]
  
  # Compute the pairwise correlation matrix
  cor_matrix <- cor(t(dataCPM_sub), method = "pearson")
  mapping <- read.csv("Gene_name and Protein_name.csv", stringsAsFactors = FALSE)
  # Remove extra spaces and standardize names to lowercase for matching
  mapping[, 1] <- tolower(trimws(mapping[, 1]))
  gene_names <- tolower(trimws(gene_names))
  gene_names <- mapping[match(gene_names, mapping[, 1]), 2]
  rownames(cor_matrix) <- gene_names
  colnames(cor_matrix) <- gene_names
  
  # Generate the correlation plot and save as PNG
  file_name <- paste0("corrplot_", project, ".png")
  png(filename = file_name, width = 1000, height = 800)
  corrplot(cor_matrix, method = "color", type = "lower",  title = paste("Correlation Comparison for", project),  mar = c(0, 0, 1, 0), 
           tl.col = "black", tl.srt = 45)
  dev.off()
  
  # Generate the correlation plot and save as PNG
  file_name <- paste0("corrplot_", project, "_AOE", ".png")
  png(filename = file_name, width = 1000, height = 800)
  corrplot(cor_matrix, method = "color", type = "lower", order = "AOE", title = paste("Correlation Comparison for", project),  mar = c(0, 0, 1, 0),
           tl.col = "black", tl.srt = 45)
  dev.off()
}