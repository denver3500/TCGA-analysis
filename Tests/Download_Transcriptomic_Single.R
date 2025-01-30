################################# Libraries
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(gridExtra)
library(ggsignif)

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

# Retrieve project information
projects_info <- TCGAbiolinks:::getGDCprojects()
project_id <- "TCGA-CHOL"
project_name <- projects_info$name[projects_info$project_id == project_id]
disease_type <- projects_info$disease_type[projects_info$project_id == project_id]

plot_list <- list() # List to store plots

for (gene in rownames(filtered_TPM_data)) { # Cycle for each gene
  gene_data <- data.frame(
    sample = colnames(filtered_TPM_data),
    expression = filtered_TPM_data[gene, ],
    type = sample_types$type
  )
  
  readable_name <- genes_filter$Gene.name[genes_filter$ENSEMBL == gene] # Substitute ENSEMBL id with Gene name
  
  # Subset data for "Solid Tissue Normal" and "Primary Tumor"
  normal_data <- gene_data$expression[gene_data$type == "Solid Tissue Normal"]
  tumor_data <- gene_data$expression[gene_data$type == "Primary Tumor"]
  
  # Perform Wilcoxon test
  wilcox_test <- wilcox.test(tumor_data, normal_data)
  
  # Determine y-axis limits with padding
  max_expression <- max(gene_data$expression)
  y_limit <- max_expression + (max_expression * 0.1) # Add 10% padding
  
  # Create the plot
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
    ) +
    ylim(0, y_limit) # Set y-axis limits
  
  # Add significance annotation if the p-value is significant
  if (wilcox_test$p.value < 0.05) {
    p <- p + geom_signif(
      comparisons = list(c("Solid Tissue Normal", "Primary Tumor")),
      map_signif_level = TRUE
    )
  }
  
  plot_list[[readable_name]] <- p # Write plot to plot_list
}

# Calculate font size based on the length of the text
text_length <- nchar(paste("Disease:", disease_type))
font_size <- ifelse(text_length > 50, 4, 6) # Adjust font size based on text length

# Wrap the text to fit within the plot
wrapped_text <- paste(strwrap(paste("Disease:", disease_type), width = 30), collapse = "\n")

# Create a plot for the project description
description_plot <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = wrapped_text, 
           hjust = 0.5, vjust = 0.5, size = font_size, color = "black") +
  theme_void() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# Create an annotation plot for the Wilcoxon test information
annotation_plot <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Wilcoxon test\n* p < 0.05\n** p < 0.01\n*** p < 0.001", 
           hjust = 0.5, vjust = 0.5, size = 5, color = "black") +
  theme_void() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# Create a plot for the project name
project_name_plot <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = project_name, 
           hjust = 0.5, vjust = 0.5, size = 6, color = "black") +
  theme_void() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# Calculate the number of rows needed
num_plots <- length(plot_list)
num_cols <- 4
num_rows <- ceiling((num_plots + 3) / num_cols) # +3 for the project name, description, and annotation plots

# Create the layout matrix
layout_matrix <- matrix(NA, nrow = num_rows, ncol = num_cols)
layout_matrix[1:num_plots] <- 1:num_plots
layout_matrix[num_plots + 1] <- num_plots + 1 # Project name plot
layout_matrix[num_plots + 2] <- num_plots + 2 # Description plot
layout_matrix[num_plots + 3] <- num_plots + 3 # Annotation plot

# Combine all plots into one image
combined_plot <- grid.arrange(
  grobs = c(plot_list, list(project_name_plot, description_plot, annotation_plot)),
  ncol = num_cols,
  layout_matrix = layout_matrix
)

# Save the combined plot
ggsave(filename = "CHOL_genes_plot.png", plot = combined_plot, width = 20, height = 16)