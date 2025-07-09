library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(gridExtra)
library(ggsignif)

rds_files <- list.files("raw_data/", pattern = "*.rds") # Looking for created .rds files after downloading data
filtered_gene_data <- list() # Store data for later
primary_normal_data <- list() # Store primary tumor and solid tissue normal data for all projects # nolint: line_length_linter.
genes_filter <- read.csv("TCGA-gene_expression/ProteinID.csv") # File with Gene name, UniProt ID, AGID (For Proteomics), and ENSEMBL id
Ensembl_column <- genes_filter$ENSEMBL
projects_info <- TCGAbiolinks:::getGDCprojects() # Load information about projects

for (rds_file in rds_files) { 
  message(paste("Processing file:", rds_file))
  transcriptomic_exp <- readRDS(file.path("raw_data", rds_file))
  
  if (!inherits(transcriptomic_exp, "SummarizedExperiment")) { # Check if the loaded object is a SummarizedExperiment
    message(paste("Error: The file", rds_file, "does not contain a SummarizedExperiment object."))
    next
  }
  
  project_id <- gsub("_transcriptomic_exp.rds", "", rds_file) %>% gsub("raw_data/", "", .) # Extract project ID from the file name
  project_name <- projects_info$name[projects_info$project_id == project_id] #Retrieve project name
  disease_type <- projects_info$disease_type[projects_info$project_id == project_id] #Retrieve disease type
  sanitized_project_name <- gsub("[^A-Za-z0-9]", "_", project_name) # Make sure that project name .png can be saved
  
  TPM_data <- assay(transcriptomic_exp, 'tpm_unstrand')
  TPM_data <- TPM_data[rownames(TPM_data) %in% genes_filter$ENSEMBL, ]
  filtered_gene_data[[project_id]] <- TPM_data
  sample_types <- colData(transcriptomic_exp) %>% as.data.frame()
  sample_types <- sample_types[colnames(TPM_data), ] # Subset sample_types to include only the samples present in filtered_TPM_data
  
  if ("sample_type" %in% colnames(sample_types)) {
    sample_types$type <- factor(sample_types$sample_type, levels = c("Solid Tissue Normal", "Primary Tumor")) # Change order so Solid Tissue Normal is first
  } else {
    message(paste("Error: 'sample_type' column not found in project:", project_id))
    next
  }
  
  project_primary_normal_data <- list()   # Create a list to store primary tumor and solid tissue normal data for the current project
  plot_list <- list()
  
  for (gene in rownames(TPM_data)) { # Cycle for each gene
    gene_data <- data.frame(
      sample = colnames(TPM_data),
      expression = TPM_data[gene, ],
      type = sample_types$type
    )
    
    readable_name <- genes_filter$Gene.name[genes_filter$ENSEMBL == gene] # Substitute ENSEMBL id with Gene name
    
    # Subset data for "Solid Tissue Normal" and "Primary Tumor"
    normal_data <- gene_data$expression[gene_data$type == "Solid Tissue Normal"]
    tumor_data <- gene_data$expression[gene_data$type == "Primary Tumor"]
    
    # Store the primary tumor and solid tissue normal data in the list
    project_primary_normal_data[[readable_name]] <- list(
      normal_data = normal_data,
      tumor_data = tumor_data
    )
    
    # Check if both groups have enough observations
    if (length(normal_data) > 1 && length(tumor_data) > 1) {
      # Perform Wilcoxon test
      wilcox_test <- wilcox.test(tumor_data, normal_data)
      
      # Determine y-axis limits with padding
      max_expression <- max(gene_data$expression)
      y_limit <- max_expression + (max_expression * 0.1) # Add 10% padding
      
      # Create the plot
      p <- ggplot(gene_data, aes(x = type, y = expression, fill = type)) +
        geom_boxplot() +
        labs(title = paste("Expression Levels of", readable_name, "in", project_id),
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
      
      plot_list[[paste(project_id, readable_name, sep = "_")]] <- p # Write plot to plot_list
    } else {
      message(paste("Not enough observations for Wilcoxon test in project:", project_id, "gene:", readable_name))
    }
  }
  
  # Store the primary tumor and solid tissue normal data for the current project
  primary_normal_data[[project_id]] <- project_primary_normal_data
  
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
  
  # Check if plot_list is not empty before calling grid.arrange
  if (length(plot_list) > 0) {
    # Calculate the number of rows needed
    num_plots <- length(plot_list)
    num_cols <- 4
    num_rows <- ceiling((num_plots + 2) / num_cols) # +2 for the project name and annotation plots
    
    # Create the layout matrix
    layout_matrix <- matrix(NA, nrow = num_rows, ncol = num_cols)
    layout_matrix[1:num_plots] <- 1:num_plots
    layout_matrix[num_plots + 1] <- num_plots + 1 # Project name plot
    layout_matrix[num_plots + 2] <- num_plots + 2 # Annotation plot
    
    # Combine all plots for the current project into one image
    combined_plot <- grid.arrange(
      grobs = c(plot_list, list(project_name_plot, annotation_plot)),
      ncol = num_cols,
      layout_matrix = layout_matrix
    )
    
    # Calculate dynamic width and height
    plot_width <- num_cols * 5 # 5 inches per column
    plot_height <- num_rows * 5 # 5 inches per row
    
    # Save the combined plot for the current project using the sanitized project name
    ggsave(
      filename = file.path("pictures", "expression_projects", paste0(project_id, "_", sanitized_project_name, "_genes_plot.png")),
      plot = combined_plot,
      width = plot_width,
      height = plot_height
    )
  } else {
    message(paste("No valid plots for project:", project_id))
  }
}

# Save the filtered gene data for all projects
saveRDS(filtered_gene_data, file = "TCGA-gene_expression/filtered_gene_data.rds")

# Save the primary tumor and solid tissue normal data for all projects
saveRDS(primary_normal_data, file = "TCGA-gene_expression/primary_normal_data.rds")
