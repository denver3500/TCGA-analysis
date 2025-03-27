library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(gridExtra)
library(ggsignif)

rds_files <- list.files("TCGA-Chaperones/rds", pattern = "*.rds", full.names = TRUE) # Looking for created .rds files after downloading data
filtered_gene_data <- list() # Store data for later
primary_normal_data <- list() # Store primary tumor and solid tissue normal data for all projects # nolint: line_length_linter.
genes_filter <- read.csv("TCGA-Chaperones/gene_list.csv") # File with Gene name, UniProt ID, AGID (For Proteomics), and ENSEMBL id
Ensembl_column <- genes_filter$ENSEMBL
projects_info <- TCGAbiolinks:::getGDCprojects() # Load information about projects

for (rds_file in rds_files) {
  message(paste("Processing file:", rds_file))
  transcriptomic_exp <- readRDS(rds_file) 

  if (!inherits(transcriptomic_exp, "SummarizedExperiment")) { # Check if the loaded object is a SummarizedExperiment
    message(paste("Error: The file", rds_file, "does not contain a SummarizedExperiment object."))
    next
  }

project_id <- gsub(".*/(.*?)_.*\\.rds$", "\\1", rds_file)  # Extracts just the project ID from the full path
  project_name <- projects_info$name[projects_info$project_id == project_id] #Retrieve project name
  disease_type <- projects_info$disease_type[projects_info$project_id == project_id] #Retrieve disease type
  sanitized_project_name <- gsub("[^A-Za-z0-9]", "_", project_name) # Make sure that project name .png can be saved

  TPM_data <- assay(transcriptomic_exp, "tpm_unstrand")
  TPM_data <- TPM_data[sapply(rownames(TPM_data), function(x) any(sapply(genes_filter$ENSEMBL, function(pattern) grepl(pattern, x)))), ]
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
  # Clear plot list for this project
  plot_list <- list()
  
  for (gene in rownames(TPM_data)) { # Cycle for each gene
    gene_data <- data.frame(
      sample = colnames(TPM_data),
      expression = TPM_data[gene, ],
      type = sample_types$type
    )
    
    matching_idx <- which(sapply(genes_filter$ENSEMBL, function(pattern) grepl(pattern, gene)))[1]
    if (!is.na(matching_idx) && length(matching_idx) > 0) {
      readable_name <- genes_filter$Name[matching_idx]
    } else {
      readable_name <- gene  # If no match is found, use the original ID
    }
    
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
      max_expression <- max(gene_data$expression, na.rm = TRUE)
      y_limit <- max_expression * 1.2  # Add 20% padding for significance bars
      
      # Get significance level for annotation
      sig_level <- if(wilcox_test$p.value < 0.001) {
        "***"
      } else if(wilcox_test$p.value < 0.01) {
        "**"
      } else if(wilcox_test$p.value < 0.05) {
        "*"
      } else {
        "ns"
      }
      
      # Create the plot with larger dimensions
      p <- ggplot(gene_data, aes(x = type, y = expression, fill = type)) +
        geom_boxplot(width = 0.6, outlier.size = 1.5) +
        scale_fill_manual(values = c("Solid Tissue Normal" = "#619CFF", "Primary Tumor" = "#F8766D")) +
        labs(title = paste0(readable_name, " (", sig_level, ")"),
             x = NULL,
             y = "TPM") +
        theme_minimal() +
        theme(
          panel.background = element_rect(fill = "white", colour = "grey90"),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          plot.margin = unit(c(0.5, 0.5, 1.5, 0.5), "cm")
        ) +
        ylim(0, y_limit)
      
      # Add significance annotation if the p-value is significant
      if (wilcox_test$p.value < 0.05) {
        p <- p + geom_signif(
          comparisons = list(c("Solid Tissue Normal", "Primary Tumor")),
          map_signif_level = TRUE,
          tip_length = 0.02,
          vjust = 0.5
        )
      }
      
      plot_list[[readable_name]] <- p
    } else {
      message(paste("Not enough observations for Wilcoxon test in project:", project_id, "gene:", readable_name))
    }
  }
  
  # Store the primary tumor and solid tissue normal data for the current project
  primary_normal_data[[project_id]] <- project_primary_normal_data
  
  # Create a project name plot with proper text wrapping
  create_project_name_plot <- function(project_name) {
    # Use stringr to wrap text
    if(requireNamespace("stringr", quietly = TRUE)) {
      if(nchar(project_name) > 30) {
        project_name <- stringr::str_wrap(project_name, width = 30)
      }
    }
    
    # Determine font size based on name length
    font_size <- 16
    if(nchar(project_name) > 50) font_size <- 12
    if(nchar(project_name) > 80) font_size <- 10
    
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = project_name, 
               hjust = 0.5, vjust = 0.5, size = font_size, face = "bold", color = "black") +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
      )
  }
  
  # Create the project name plot
  project_name_plot <- create_project_name_plot(project_name)
  
  # Create a Wilcoxon test legend
  wilcox_legend <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "Wilcoxon rank sum test\n* p < 0.05\n** p < 0.01\n*** p < 0.001\nns: not significant", 
             hjust = 0.5, vjust = 0.5, size = 10, color = "black") +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  # Add sample size information
  sample_info <- paste0(
    "Sample sizes:\n",
    "Normal tissue: ", sum(sample_types$type == "Solid Tissue Normal"), "\n",
    "Tumor tissue: ", sum(sample_types$type == "Primary Tumor")
  )
  
  sample_info_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = sample_info, 
             hjust = 0.5, vjust = 0.5, size = 12, color = "black") +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  # Check if plot_list is not empty
  if (length(plot_list) > 0) {
    # Calculate layout dimensions based on number of plots
    num_plots <- length(plot_list)
    
    # Determine optimal grid layout
    if (num_plots <= 4) {
      num_cols <- 2
    } else if (num_plots <= 9) {
      num_cols <- 3
    } else {
      num_cols <- 4
    }
    
    num_rows <- ceiling(num_plots / num_cols)
    
    # Create a layout that includes project name and legend
    # Top row for project name, last row for Wilcoxon legend and sample info
    layout_matrix <- matrix(NA, nrow = num_rows + 2, ncol = num_cols)
    
    # First row is for project name, spanning all columns
    layout_matrix[1, ] <- rep(num_plots + 1, num_cols)
    
    # Fill in the gene plots
    for (i in 1:num_plots) {
      row_idx <- floor((i-1) / num_cols) + 2  # +2 because first row is project name
      col_idx <- ((i-1) %% num_cols) + 1
      layout_matrix[row_idx, col_idx] <- i
    }
    
    # Last row, first half for Wilcoxon legend
    layout_matrix[num_rows + 2, 1:floor(num_cols/2)] <- rep(num_plots + 2, floor(num_cols/2))
    
    # Last row, second half for sample info
    layout_matrix[num_rows + 2, (floor(num_cols/2) + 1):num_cols] <- rep(num_plots + 3, ceiling(num_cols/2))
    
    # Get list of gene plots in the right order
    ordered_plots <- plot_list[names(plot_list)]
    
    # Combine with title and legend plots
    all_plots <- c(ordered_plots, 
                  list(project_name_plot, wilcox_legend, sample_info_plot))
    
    # Create the complete plot
    dir.create("TCGA-Chaperones/pictures/expression_projects", showWarnings = FALSE, recursive = TRUE)
    
    # Calculate dimensions based on number of plots
    plot_width <- num_cols * 5  # 5 inches per column
    plot_height <- (num_rows + 2) * 4  # 4 inches per row
    
    # Use gridExtra to arrange the plots
    combined_plot <- gridExtra::grid.arrange(
      grobs = all_plots,
      layout_matrix = layout_matrix,
      widths = rep(1, num_cols),
      heights = c(0.8, rep(1, num_rows), 0.8)  # Smaller height for title and legend rows
    )
    
    # Save the plot
    ggsave(
      filename = file.path("TCGA-Chaperones/pictures/expression_projects", 
                      paste0(project_id, "_", sanitized_project_name, "_gene_expression.png")),
      plot = combined_plot,
      width = plot_width,
      height = plot_height,
      dpi = 300
    )

    message(paste("Created plot for project:", project_id, "with", num_plots, "genes"))
  } else {
    message(paste("No valid plots for project:", project_id))
  }
}

# Save the filtered gene data for all projects
saveRDS(filtered_gene_data, file = "TCGA-Chaperones/filtered_gene_data.rds")

# Save the primary tumor and solid tissue normal data for all projects
saveRDS(primary_normal_data, file = "TCGA-Chaperones/primary_normal_data.rds")
higher_in_tumor_counts <- data.frame(
  project_id = character(),
  higher_in_tumor_count = integer(),
  stringsAsFactors = FALSE
)

# Process each project's data
for (project_id in names(primary_normal_data)) {
  project_data <- primary_normal_data[[project_id]]
  higher_in_tumor_count <- 0
  total_genes_compared <- 0
  
  for (gene_name in names(project_data)) {
    gene_data <- project_data[[gene_name]]
    
    # Check if we have sufficient data for both normal and tumor
    if (length(gene_data$normal_data) > 1 && length(gene_data$tumor_data) > 1) {
      total_genes_compared <- total_genes_compared + 1
      
      # Calculate median expression for normal and tumor
      median_normal <- median(gene_data$normal_data, na.rm = TRUE)
      median_tumor <- median(gene_data$tumor_data, na.rm = TRUE)
      
      # Perform Wilcoxon test to check for statistical significance
      wilcox_result <- wilcox.test(gene_data$tumor_data, gene_data$normal_data)
      p_value <- wilcox_result$p.value
      
      # Only count if expression is higher in tumor AND statistically significant
      if (median_tumor > median_normal && p_value < 0.05) {
        higher_in_tumor_count <- higher_in_tumor_count + 1
      }
    }
  }
  
  # Only add projects that have at least one gene with significantly higher expression in tumor
  if (higher_in_tumor_count > 0) {
    # Add to results dataframe
    new_row <- data.frame(
      project_id = project_id,
      higher_in_tumor_count = higher_in_tumor_count,
      stringsAsFactors = FALSE
    )
    
    higher_in_tumor_counts <- rbind(higher_in_tumor_counts, new_row)
  }
}

# Write the results to a CSV file
write.csv(higher_in_tumor_counts, 
          file = "TCGA-Chaperones/significant_higher_in_tumor_gene_counts.csv", 
          row.names = FALSE)

# Print a summary to console
message(paste("Analysis complete. Results saved to TCGA-Chaperones/significant_higher_in_tumor_gene_counts.csv"))
message(paste("Total projects with genes significantly higher in tumor:", nrow(higher_in_tumor_counts)))