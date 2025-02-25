################################# Libraries
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(gridExtra)
library(ggsignif)

primary_normal_data <- readRDS(file.path("raw_data", "primary_normal_data.rds")) # Load the file we got running Process transcriptomic.R
projects_info <- TCGAbiolinks:::getGDCprojects() # Load information about projects
project_id_to_name <- setNames(projects_info$name, projects_info$project_id) # Create a mapping between project IDs and project names

# List of genes of interest
genes_of_interest <- c("Hsp90β", "p23", "FKBP52", "Hsp90⍺", "FKBP51", "GR", "CDK4", "Hop", "Hsp70")

# Create a data frame to store the top 10 most expressed projects for each gene
if (!exists("top_expressed_projects")) {
  top_expressed_projects <- data.frame()
}

for (gene_of_interest in genes_of_interest) {
  primary_tumor_data <- data.frame()
  
  for (project_id in names(primary_normal_data)) {
    project_data <- primary_normal_data[[project_id]]
    
    if (gene_of_interest %in% names(project_data)) {
      gene_data <- project_data[[gene_of_interest]]
      tumor_data <- gene_data$tumor_data
      
      # Combine the data into a single data frame
      primary_tumor_data <- rbind(
        primary_tumor_data,
        data.frame(
          project = project_id,
          expression = tumor_data
        )
      )
    }
  }
  
  primary_tumor_data$project <- project_id_to_name[primary_tumor_data$project] # Substitute project IDs with project names
  
  # Sort projects by median expression
  primary_tumor_data <- primary_tumor_data %>%
    group_by(project) %>%
    mutate(median_expression = median(expression)) %>%
    ungroup() %>%
    arrange(median_expression) %>%
    mutate(project = factor(project, levels = unique(project)))
  
  # Collect top 10 most expressed projects for this gene
  if (nrow(primary_tumor_data) > 0) {
    top10 <- primary_tumor_data %>%
      group_by(project) %>%
      summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
      arrange(desc(mean_expression)) %>%
      slice_head(n = 10) %>%
      mutate(gene = gene_of_interest)
    
    # Append to a global data frame
    top_expressed_projects <- rbind(top_expressed_projects, top10)
    
    # Calculate y-axis limits
    y_limits <- range(primary_tumor_data$expression, na.rm = TRUE)
    
    # Calculate plot dimensions
    num_projects <- length(unique(primary_tumor_data$project))
    max_project_name_length <- max(nchar(as.character(primary_tumor_data$project)))
    plot_width <- max(10, num_projects * 0.5) + max_project_name_length * 0.2
    plot_height <- 8 + max_project_name_length * 0.1
    
    # Create the plot
    primary_tumor_plot <- ggplot(
      primary_tumor_data,
      aes(x = project, y = expression, fill = project)
    ) +
      geom_boxplot() +
      labs(
        title = paste("Primary Tumor Expression Levels of", gene_of_interest),
        x = "Project",
        y = "Expression (TPM)"
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(6, 6, 6, 6), "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)
      ) +
      ylim(y_limits)
    
    # Display the plot
    print(primary_tumor_plot)
    
    # Save the plot
    ggsave(
      filename = file.path(
        "pictures",
        "expression_amongsts_projects",
        paste0("Primary_Tumor_Expression_Levels_", gene_of_interest, ".png")
      ),
      plot = primary_tumor_plot,
      width = plot_width,
      height = plot_height
    )
  } else {
    message(paste("No data available for gene:", gene_of_interest))
  }
}

# Create the table from the project column in top_expressed_projects
unique_counts <- table(as.vector(top_expressed_projects$project))
print(unique_counts)

# Convert the table object to a data frame for easier CSV export
unique_counts_df <- as.data.frame(unique_counts)

# Sort by frequency (the "Freq" column) in descending order
unique_counts_df <- unique_counts_df[order(-unique_counts_df$Freq), ]

# Write the data frame to a CSV file in the "analysis" folder
write.csv(unique_counts_df, file.path("analysis", "unique_counts.csv"), row.names = FALSE)
print(unique_counts_df)