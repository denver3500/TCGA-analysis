library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)

# Read the gene counts and filter projects with more than 10 significantly increased genes
gene_counts <- read.csv("TCGA-Chaperones/significant_higher_in_tumor_gene_counts.csv", stringsAsFactors = FALSE)
filtered_projects <- gene_counts %>% 
  filter(higher_in_tumor_count >= 10) %>%
  pull(project_id)

message(paste("Found", length(filtered_projects), "projects with 10+ significantly increased genes"))

# Load the primary_normal_data which contains the expression values
primary_normal_data <- readRDS("TCGA-Chaperones/primary_normal_data.rds")

# Filter the data to only include selected projects
primary_normal_data <- primary_normal_data[filtered_projects]

# Get project names from TCGAbiolinks for better labels
projects_info <- TCGAbiolinks:::getGDCprojects()
project_labels <- sapply(filtered_projects, function(project) {
  name <- projects_info$name[projects_info$project_id == project]
  if(length(name) == 0) return(project) # Use project ID if name not found
  return(paste0(project, ": ", name))
})

# Get the list of all genes
all_genes <- unique(unlist(lapply(primary_normal_data, names)))
message(paste("Found", length(all_genes), "unique genes across all projects"))

# Create directory for the output
dir.create("TCGA-Chaperones/pictures/cross_project_expression", showWarnings = FALSE, recursive = TRUE)

# For each gene, create a plot of expression across projects
for(gene in all_genes) {
  # Collect median tumor and normal expression for this gene across all projects
  gene_data <- data.frame(
    project = character(),
    tissue_type = character(),
    median_expression = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(project in filtered_projects) {
    if(gene %in% names(primary_normal_data[[project]])) {
      # Get tumor data
      tumor_data <- primary_normal_data[[project]][[gene]]$tumor_data
      if(length(tumor_data) > 0) {
        gene_data <- rbind(gene_data, data.frame(
          project = project,
          tissue_type = "Tumor",
          median_expression = median(tumor_data, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
      
      # Get normal data
      normal_data <- primary_normal_data[[project]][[gene]]$normal_data
      if(length(normal_data) > 0) {
        gene_data <- rbind(gene_data, data.frame(
          project = project,
          tissue_type = "Normal",
          median_expression = median(normal_data, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Skip if not enough data
  if(nrow(gene_data) < 3) {
    message(paste("Skipping", gene, "- insufficient data"))
    next
  }
  
  # Get the order of projects based on tumor expression (ascending)
  project_order <- gene_data %>%
    filter(tissue_type == "Tumor") %>%
    arrange(median_expression) %>%
    pull(project)
  
  # Add project names and preserve the order based on tumor expression
  gene_data <- gene_data %>%
    mutate(project_name = sapply(project, function(p) {
      idx <- which(filtered_projects == p)
      if(length(idx) > 0) return(project_labels[idx])
      return(p)
    })) %>%
    mutate(project_name = factor(project_name, 
                                 levels = sapply(project_order, function(p) {
                                   idx <- which(filtered_projects == p)
                                   if(length(idx) > 0) return(project_labels[idx])
                                   return(p)
                                 })))
  
  # Create the plot
  p <- ggplot(gene_data, aes(x = project_name, y = median_expression, fill = tissue_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    coord_flip() + # Horizontal bars
    scale_fill_manual(values = c("Normal" = "#619CFF", "Tumor" = "#F8766D"),
                      name = "Tissue Type") +
    labs(
      title = paste("Median Expression of", gene, "Across Cancer Projects"),
      x = "Cancer Project",
      y = "Median Gene Expression (TPM)"
    ) +
    theme_minimal() +
    theme(
      # Set explicit white backgrounds
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      # Remove horizontal grid lines
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      # Other formatting
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text.y = element_text(size = 10)
    )
  
  # Save the plot
  ggsave(
    filename = paste0("TCGA-Chaperones/pictures/cross_project_expression/", gene, "_cross_project.png"),
    plot = p,
    width = 14,
    height = max(8, nrow(gene_data) * 0.25), # Scale height based on number of projects
    dpi = 300
  )
  
  message(paste("Created cross-project expression plot for", gene))
}

message("Cross-project expression analysis complete")