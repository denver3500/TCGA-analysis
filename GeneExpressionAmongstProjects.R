################################# Libraries
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(gridExtra)
library(ggsignif)

primary_normal_data <- readRDS("primary_normal_data.rds") # Load the file we got running Process transcriptomic.R
projects_info <- TCGAbiolinks:::getGDCprojects() # Load information about projects
project_id_to_name <- setNames(projects_info$name, projects_info$project_id)#Create a mapping between project IDs and project names

# List of genes of interest
genes_of_interest <- c("Hsp90β", "p23", "FKBP52", "Hsp90⍺", "FKBP51", "GR", "CDK4", "Hop", "Hsp70") # List of genes. ENSEMBL ID'S has been replaced in Process transcriptomic.R

for (gene_of_interest in genes_of_interest) {
  primary_tumor_data <- data.frame()
  
  for (project_id in names(primary_normal_data)) {
    project_data <- primary_normal_data[[project_id]]
    
    if (gene_of_interest %in% names(project_data)) {
      gene_data <- project_data[[gene_of_interest]]
      tumor_data <- gene_data$tumor_data
      
      # Combine the data into a single data frame
      primary_tumor_data <- rbind(primary_tumor_data, data.frame(
        project = project_id,
        expression = tumor_data
      ))
    }
  }

  primary_tumor_data$project <- project_id_to_name[primary_tumor_data$project]   # Substitute project IDs with project names for better readability
  

  primary_tumor_data <- primary_tumor_data %>%   # Sort from min to max. It puts Exceptional Responders last every time for no reason whatsoever. Idk, should be related to N/A sample type that I somehow downloaded?
    group_by(project) %>%
    mutate(median_expression = median(expression)) %>%
    ungroup() %>%
    arrange(median_expression) %>%
    mutate(project = factor(project, levels = unique(project)))
  
  if (nrow(primary_tumor_data) > 0) { #Not every project has a data to plot
    # Calculate y-axis limits
    y_limits <- range(primary_tumor_data$expression, na.rm = TRUE)
    
    # Calculate plot dimensions
    num_projects <- length(unique(primary_tumor_data$project))
    max_project_name_length <- max(nchar(as.character(primary_tumor_data$project)))
    plot_width <- max(10, num_projects * 0.5) + max_project_name_length * 0.2 # Adjust width based on number of projects and max project name length
    plot_height <- 8 + max_project_name_length * 0.1 # Adjust height based on max project name length
    
    # Create the plot
    primary_tumor_plot <- ggplot(primary_tumor_data, aes(x = project, y = expression, fill = project)) +
      geom_boxplot() +
      labs(title = paste("Primary Tumor Expression Levels of", gene_of_interest),
           x = "Project",
           y = "Expression (TPM)") +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24), # Increase title font size
        plot.margin = unit(c(6, 6, 6, 6), "cm"), # Increase margins around the plot
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # Rotate x-axis text labels and increase size
        axis.text.y = element_text(size = 16), # Increase y-axis text size
        axis.title.x = element_text(size = 20), # Increase x-axis title font size
        axis.title.y = element_text(size = 20)  # Increase y-axis title font size
      ) +
      ylim(y_limits) # Set y-axis limits
    
    # Display the plot
    print(primary_tumor_plot)
    
    # Save the plot with adjusted width and height
    ggsave(filename = paste0("Primary_Tumor_Expression_Levels_", gene_of_interest, ".png"), plot = primary_tumor_plot, width = plot_width, height = plot_height)
  } else {
    message(paste("No data available for gene:", gene_of_interest))
  }
}