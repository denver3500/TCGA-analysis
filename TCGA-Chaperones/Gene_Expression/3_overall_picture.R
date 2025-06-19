library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggrepel)  # For non-overlapping labels
library(patchwork) # For combining plots
library(TCGAbiolinks) # For project metadata

# Load the saved primary_normal_data
primary_normal_data <- readRDS("TCGA-Chaperones/primary_normal_data.rds")

# Load gene family information
gene_list <- read_csv("TCGA-Chaperones/gene_list.csv", show_col_types = FALSE)

# Get project information for better labeling
projects_info <- TCGAbiolinks:::getGDCprojects()

# Process the data into a tidy format for plotting
expression_summary <- tibble()

for (project_id in names(primary_normal_data)) {
  project_data <- primary_normal_data[[project_id]]
  project_name <- projects_info$name[projects_info$project_id == project_id]
  disease_type <- projects_info$disease_type[projects_info$project_id == project_id]
  
  for (gene_name in names(project_data)) {
    gene_data <- project_data[[gene_name]]
    
    # Skip if insufficient data
    if (length(gene_data$normal_data) < 3 || length(gene_data$tumor_data) < 3) {
      next
    }
    
    # Calculate statistics
    median_normal <- median(gene_data$normal_data, na.rm = TRUE)
    median_tumor <- median(gene_data$tumor_data, na.rm = TRUE)
    mean_normal <- mean(gene_data$normal_data, na.rm = TRUE)
    mean_tumor <- mean(gene_data$tumor_data, na.rm = TRUE)
    log2fc <- log2((mean_tumor + 0.1) / (mean_normal + 0.1))  # Add small value to prevent division by zero
    
    # Statistical test
    wilcox_result <- wilcox.test(gene_data$tumor_data, gene_data$normal_data)
    p_value <- wilcox_result$p.value
    
    # Calculate effect size and confidence intervals
    n_normal <- length(gene_data$normal_data)
    n_tumor <- length(gene_data$tumor_data)
    
    # Standard error of log2FC (approximate method)
    se_normal <- sd(gene_data$normal_data, na.rm = TRUE) / sqrt(n_normal)
    se_tumor <- sd(gene_data$tumor_data, na.rm = TRUE) / sqrt(n_tumor)
    se_log2fc <- sqrt((se_normal/(mean_normal+0.1))^2 + (se_tumor/(mean_tumor+0.1))^2) / log(2)
    
    # 95% confidence intervals
    ci_low <- log2fc - 1.96 * se_log2fc
    ci_high <- log2fc + 1.96 * se_log2fc
    
    # Gather all data
    project_gene_summary <- tibble(
      project_id = project_id,
      project_name = if(length(project_name) > 0) project_name else project_id,
      disease_type = if(length(disease_type) > 0) disease_type else "Unknown",
      gene_name = gene_name,
      median_normal = median_normal,
      median_tumor = median_tumor,
      mean_normal = mean_normal,
      mean_tumor = mean_tumor,
      log2fc = log2fc,
      ci_low = ci_low,
      ci_high = ci_high,
      p_value = p_value,
      n_normal = n_normal,
      n_tumor = n_tumor,
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      is_significant = p_value < 0.05
    )
    
    expression_summary <- bind_rows(expression_summary, project_gene_summary)
  }
}

# Create a directory for the output
dir.create("TCGA-Chaperones/pictures/pan_cancer_analysis", showWarnings = FALSE, recursive = TRUE)

# Filter projects with sufficient data for most genes
gene_counts <- expression_summary %>%
  group_by(project_id) %>%
  summarise(gene_count = n_distinct(gene_name))

valid_projects <- gene_counts %>%
  filter(gene_count >= 5) %>%  # At least 5 genes per project
  pull(project_id)

filtered_summary <- expression_summary %>%
  filter(project_id %in% valid_projects)

message("Projects with sufficient data: ", length(valid_projects), " out of ", length(unique(expression_summary$project_id)))

# List of all genes in dataset
all_genes <- unique(filtered_summary$gene_name)
message("Analyzing ", length(all_genes), " chaperone genes across ", length(valid_projects), " cancer types")

# Create mapping between gene names and families
gene_family_mapping <- setNames(gene_list$Family, gene_list$Name)

# Get project name mapping
project_name_mapping <- filtered_summary %>%
  select(project_id, project_name) %>%
  distinct() %>%
  deframe()

# Create a matrix of log2FC values using project names instead of IDs
# Convert to regular data.frame before setting rownames - FIX FOR ROW NAMES ISSUE
heatmap_df <- filtered_summary %>%
  select(project_id, project_name, gene_name, log2fc) %>%
  distinct() %>%
  pivot_wider(
    id_cols = c(project_id, project_name),
    names_from = gene_name,
    values_from = log2fc
  ) %>%
  as.data.frame()  # Convert tibble to data.frame

# Set rownames properly on data.frame, not tibble
project_names <- heatmap_df$project_name
heatmap_data <- heatmap_df %>%
  select(-project_id, -project_name) %>%
  as.matrix()
rownames(heatmap_data) <- project_names

# Create a matrix of p-values for annotation - FIX ROW NAMES ISSUE
pvalue_df <- filtered_summary %>%
  select(project_id, project_name, gene_name, p_value) %>%
  distinct() %>%
  pivot_wider(
    id_cols = c(project_id, project_name),
    names_from = gene_name,
    values_from = p_value
  ) %>%
  as.data.frame()  # Convert tibble to data.frame

# Set rownames properly on data.frame
pvalue_data <- pvalue_df %>%
  select(-project_id, -project_name) %>%
  as.matrix()
rownames(pvalue_data) <- project_names

# Create annotation for significance
significance_annotation <- matrix("", nrow = nrow(pvalue_data), ncol = ncol(pvalue_data))
rownames(significance_annotation) <- rownames(pvalue_data)
colnames(significance_annotation) <- colnames(pvalue_data)

for (i in 1:nrow(pvalue_data)) {
  for (j in 1:ncol(pvalue_data)) {
    if (!is.na(pvalue_data[i, j])) {
      if (pvalue_data[i, j] < 0.001) {
        significance_annotation[i, j] <- "***"
      } else if (pvalue_data[i, j] < 0.01) {
        significance_annotation[i, j] <- "**"
      } else if (pvalue_data[i, j] < 0.05) {
        significance_annotation[i, j] <- "*"
      }
    }
  }
}

# Create column (gene) annotations for family information
gene_families <- sapply(colnames(heatmap_data), function(x) {
  if(x %in% names(gene_family_mapping)) {
    return(gene_family_mapping[x])
  } else {
    return("Unknown")
  }
})

# Create as data.frame, not tibble
col_anno <- data.frame(
  Family = gene_families,
  stringsAsFactors = FALSE
)
rownames(col_anno) <- colnames(heatmap_data)

# Create colors for gene families
family_colors <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(length(unique(col_anno$Family))),
  unique(col_anno$Family)
)

anno_colors <- list(
  Family = family_colors
)

# Create enhanced heatmap with annotations - fixed with larger fonts and more spacing
pdf("TCGA-Chaperones/pictures/pan_cancer_analysis/chaperone_expression_heatmap_enhanced.pdf", width = 24, height = 20)
pheatmap(
  heatmap_data,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  display_numbers = significance_annotation,
  fontsize_number = 16,  # Increased size for significance stars
  annotation_col = col_anno,  # Gene family annotations
  annotation_colors = anno_colors,
  main = "Chaperone Expression Changes Across TCGA Cancer Types",
  fontsize = 16,  # Main title size
  fontsize_row = 16,  # Increased size for project names
  fontsize_col = 16,  # Increased size for gene names
  angle_col = 45,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  border_color = NA,
  show_rownames = TRUE,  # Explicitly show project names
  treeheight_row = 150,  # Increase space for row dendrogram
  treeheight_col = 150,  # Increase space for column dendrogram
  cellwidth = 30,        # Increased cell width for more spacing
  cellheight = 20,       # Increase cell height
  legend = TRUE,
  # Fixed consistent scale from -2 to 2
  breaks = seq(-2, 2, length.out = 101),
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2.0", "-1.0", "0.0", "+1.0", "+2.0")
)
dev.off()

# Also create a simplified version in case the complex one fails
pdf("TCGA-Chaperones/pictures/pan_cancer_analysis/chaperone_expression_heatmap_simple.pdf", width = 24, height = 20)
pheatmap(
  heatmap_data,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  fontsize = 16,
  fontsize_row = 16,
  fontsize_col = 16,
  show_rownames = TRUE,
  main = "Chaperone Expression Changes Across TCGA Cancer Types",
  # Also fix the legend in the simple version
  breaks = seq(-2, 2, length.out = 101),
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2.0", "-1.0", "0.0", "+1.0", "+2.0")
)
dev.off()

# Also increase font size in the bubble plot
bubble_plot <- filtered_summary %>%
  mutate(
    direction = case_when(
      log2fc > 0 & is_significant ~ "Up in Cancer",
      log2fc < 0 & is_significant ~ "Down in Cancer",
      TRUE ~ "Not Significant"
    )
  ) %>%
  # Use project_name instead of project_id
  ggplot(aes(x = gene_name, y = project_name, size = abs(log2fc), color = direction)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    "Up in Cancer" = "red3",
    "Down in Cancer" = "blue3",
    "Not Significant" = "gray80"
  )) +
  scale_size_continuous(range = c(0.5, 5), name = "|Log2FC|") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "right"
  ) +
  labs(
    title = "Pan-Cancer Chaperone Expression Changes",
    x = "Chaperone Gene",
    y = "Cancer Type",
    color = "Expression Change"
  )

ggsave(
  filename = "TCGA-Chaperones/pictures/pan_cancer_analysis/chaperone_bubble_plot_enhanced.pdf", 
  plot = bubble_plot,
  width = 16, 
  height = 14
)
ggsave(
  filename = "TCGA-Chaperones/pictures/pan_cancer_analysis/chaperone_bubble_plot_enhanced.png", 
  plot = bubble_plot,
  width = 16, 
  height = 14,
  dpi = 300
)

# 4. SUMMARY STATISTICS
# ---------------------
# Calculate overall statistics about chaperone expression changes

summary_stats <- filtered_summary %>%
  group_by(gene_name) %>%
  summarize(
    projects_analyzed = n(),
    up_in_cancer = sum(log2fc > 0 & is_significant),
    down_in_cancer = sum(log2fc < 0 & is_significant),
    pct_up = round(100 * up_in_cancer / projects_analyzed, 1),
    pct_down = round(100 * down_in_cancer / projects_analyzed, 1),
    mean_log2fc = mean(log2fc, na.rm = TRUE),
    median_log2fc = median(log2fc, na.rm = TRUE)
  ) %>%
  arrange(desc(pct_up))

# Add family information to summary stats
summary_stats <- summary_stats %>%
  left_join(gene_list %>% select(Name, Family), by = c("gene_name" = "Name"))

# Save the summary statistics
write_csv(summary_stats, "TCGA-Chaperones/pictures/pan_cancer_analysis/chaperone_summary_stats.csv")

# Print some summary information
message("\nChaperone Expression Summary:")
message("----------------------------")
message("Most consistently upregulated chaperones across cancer types:")
print(head(summary_stats %>% arrange(desc(pct_up)), 5))

message("\nMost consistently downregulated chaperones across cancer types:")
print(head(summary_stats %>% arrange(desc(pct_down)), 5))

message("\nAnalysis complete. Results saved to TCGA-Chaperones/pictures/pan_cancer_analysis/")