################################# Libraries
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(gridExtra)
################################# Download Transciptomic data

query.transcriptomic.chol <- GDCquery(
  project = "TCGA-CHOL", #Adenomas and Adenocarcinomas
  data.category = "Transcriptome Profiling", #Transciptomic data
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)

GDCdownload(
  query = query.transcriptomic.chol,
  files.per.chunk = 100
)
################################# Generate and save Summarizedexperiment
chol.transcriptomic.exp <- GDCprepare(
  query = query.transcriptomic.chol, 
  save = TRUE, 
  save.filename = "choltranscriptomicexp.rda"
)

################################# Filter genes of interest

genes_filter <- read.csv("ProteinID.csv") #File with Gene name, Uniprot ID, AGID (For Proteomics), and ENSEMBL id
Ensembl_column <- genes_filter$ENSEMBL
TPM_data <- assay(chol.transcriptomic.exp, 'tpm_unstrand') #Access TPM data
TPM_data_df <- as.data.frame(TPM_data)

TPM_data_df_col <- TPM_data_df %>% #Create a column with ENSEMBL ids to filter
  rownames_to_column(var = "ENSEMBL")

filtered_TPM_data_df <- TPM_data_df_col %>% #Filter by ENSEMBL id
  filter(ENSEMBL %in% Ensembl_column)

filtered_TPM_data <- filtered_TPM_data_df %>% #Remove ENSEMBL column and return to matrix form
  column_to_rownames(var = "ENSEMBL") %>%
  as.matrix()

################################# Plotting
samples.primary.tumour <- chol.transcriptomic.exp$barcode[chol.transcriptomic.exp$shortLetterCode == "TP"] #"Primary Tumor samples", Rest is "Solid Tissue Normal" since this is the only 2 sample types I've downloaded

sample_types <- data.frame( #Create a dataframe with sample type
  sample = colnames(filtered_TPM_data),
  type = ifelse(colnames(filtered_TPM_data) %in% samples.primary.tumour, "Primary Tumor", "Solid Tissue Normal")
)

sample_types$type <- factor(sample_types$type, levels = c("Solid Tissue Normal", "Primary Tumor")) #Change order so Solid Tissue Normal is first


plot_list <- list() #List to store plots

for (gene in rownames(filtered_TPM_data)) { #Cycle for each gene
  gene_data <- data.frame(
    sample = colnames(filtered_TPM_data),
    expression = filtered_TPM_data[gene, ],
    type = sample_types$type
  )
  
  readable_name <- genes_filter$Gene.name[genes_filter$ENSEMBL == gene] #Substitute ENSEMBL id with Gene name
  
  p <- ggplot(gene_data, aes(x = type, y = expression, fill = type)) +
    geom_boxplot() +
    labs(title = paste("Expression Levels of", readable_name),
         x = NULL,
         y = "TPM") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),  
      plot.margin = unit(c(1, 1, 1, 1), "cm")  
    )
  

  plot_list[[readable_name]] <- p #Write plot to plot_list
}

combined_plot <- grid.arrange( #Combine all plots into 1 image
  grobs = plot_list,
  ncol = 3,
)

ggsave(filename = "CHOL_genes_plot.png", plot = combined_plot, width = 16, height = 12) #Profit