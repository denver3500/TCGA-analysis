# TCGA Chaperone Expression-Based Survival Analysis
# This script loads all TCGA projects, performs DESeq2 normalization, 
# and creates survival plots based on chaperone gene expression

target_dir <- "TCGA-Chaperones/Survival_plots"
if (!endsWith(getwd(), target_dir)) {
  setwd(target_dir)
}

# Load required libraries
library(TCGAbiolinks)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(stringr)
library(DESeq2)
library(SummarizedExperiment)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tibble)

# Set up logging
log_file <- "chaperone_survival_analysis_log.txt"
con <- file(log_file, "w")
sink(con, split = TRUE)

message("=== TCGA Chaperone Expression-Based Survival Analysis ===")
message("Starting at: ", Sys.time())

# Create output directories
output_dir <- "pictures"
individual_dir <- file.path(output_dir, "individual_chaperones")
combined_dir <- file.path(output_dir, "combined_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(individual_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

# Load gene list
gene_list <- read_csv("../gene_list.csv", show_col_types = FALSE)
chaperone_genes <- gene_list$Name
message("Loaded ", length(chaperone_genes), " chaperone genes")

# Get all TCGA .rds files from raw_data directory
raw_data_dir <- "../../raw_data"
tcga_files <- list.files(raw_data_dir, pattern = "^TCGA-.*_transcriptomic_exp\\.rds$", full.names = TRUE)
message("Found ", length(tcga_files), " TCGA files")

# Function to extract project ID from filename
extract_project_id <- function(filename) {
  basename(filename) %>%
    str_replace("_transcriptomic_exp\\.rds$", "")
}

# Function to process a single TCGA dataset and perform DESeq2 normalization
process_tcga_for_survival <- function(file_path) {
  project_id <- extract_project_id(file_path)
  message("Processing ", project_id, "...")
  
  tryCatch({
    # Load the data
    data <- readRDS(file_path)
    
    # Check if it's a SummarizedExperiment object
    if (!inherits(data, "SummarizedExperiment")) {
      message("  Unknown data structure for ", project_id, ". Skipping.")
      return(NULL)
    }
    
    # Extract expression matrix and sample info
    expr_matrix <- assay(data)
    sample_info <- colData(data)
    gene_info <- rowData(data)
    
    # Filter for chaperone genes
    gene_matches <- character(0)
    
    # Try matching with base ENSEMBL IDs
    ensembl_ids <- str_replace(rownames(expr_matrix), "\\.\\d+$", "")
    gene_matches <- rownames(expr_matrix)[ensembl_ids %in% gene_list$ENSEMBL]
    
    # If no matches, try gene symbols
    if (length(gene_matches) == 0 && "gene_name" %in% colnames(gene_info)) {
      gene_matches <- rownames(expr_matrix)[gene_info$gene_name %in% gene_list$Name]
    }
    
    if (length(gene_matches) == 0) {
      message("  No chaperone genes found in ", project_id, ". Skipping.")
      return(NULL)
    }
    
    message("  Found ", length(gene_matches), " chaperone genes in ", project_id)
    
    # Filter expression data for chaperone genes
    chaperone_expr <- expr_matrix[gene_matches, , drop = FALSE]
    
    # Create mapping from ENSEMBL to gene names for labeling
    gene_names <- character(length(gene_matches))
    names(gene_names) <- gene_matches
    
    for (i in seq_along(gene_matches)) {
      gene_id <- gene_matches[i]
      base_ensembl <- str_replace(gene_id, "\\.\\d+$", "")
      
      # Find corresponding gene name
      if ("gene_name" %in% colnames(gene_info)) {
        gene_name <- gene_info$gene_name[rownames(expr_matrix) == gene_id]
        if (length(gene_name) > 0 && !is.na(gene_name[1])) {
          gene_names[gene_id] <- gene_name[1]
        }
      }
      
      # If no gene name found, use gene list
      if (gene_names[gene_id] == "" || is.na(gene_names[gene_id])) {
        name_match <- gene_list$Name[gene_list$ENSEMBL == base_ensembl]
        if (length(name_match) > 0) {
          gene_names[gene_id] <- name_match[1]
        } else {
          gene_names[gene_id] <- base_ensembl
        }
      }
    }
    
    # Prepare data for DESeq2 normalization
    message("  Performing DESeq2 normalization...")
    
    # Convert to integer matrix and remove any zero-only genes
    count_matrix <- round(chaperone_expr)
    count_matrix <- count_matrix[rowSums(count_matrix) > 0, , drop = FALSE]
    
    # Create sample metadata
    sample_metadata <- data.frame(
      sample_id = colnames(count_matrix),
      condition = rep("tumor", ncol(count_matrix)),
      row.names = colnames(count_matrix)
    )
    
    # Create DESeq2 object with intercept-only design
    dds <- DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = sample_metadata,
      design = ~ 1
    )
    
    # Filter low count genes
    keep <- rowSums(counts(dds) >= 10) >= ceiling(0.1 * ncol(count_matrix))
    dds <- dds[keep, ]
    
    # Estimate size factors and dispersions
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    
    # Get normalized counts
    normalized_counts <- counts(dds, normalized = TRUE)
    
    # Log2 transform (add pseudocount)
    log2_normalized <- log2(normalized_counts + 1)
    
    message("  DESeq2 normalization complete. ", nrow(log2_normalized), " genes retained.")
    
    # Create result list
    result <- list(
      project_id = project_id,
      normalized_expression = log2_normalized,
      gene_names = gene_names[rownames(log2_normalized)],
      sample_info = sample_info
    )
    
    return(result)
    
  }, error = function(e) {
    message("  Error processing ", project_id, ": ", e$message)
    return(NULL)
  })
}

# Function to create expression-based patient groups
create_expression_groups <- function(expression_vector, gene_name, method = "median") {
  if (method == "median") {
    threshold <- median(expression_vector, na.rm = TRUE)
    groups <- ifelse(expression_vector >= threshold, "High", "Low")
  } else if (method == "tertiles") {
    quantiles <- quantile(expression_vector, probs = c(1/3, 2/3), na.rm = TRUE)
    groups <- cut(expression_vector, 
                  breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
                  labels = c("Low", "Medium", "High"),
                  include.lowest = TRUE)
  } else if (method == "quartiles") {
    quantiles <- quantile(expression_vector, probs = c(0.25, 0.75), na.rm = TRUE)
    groups <- cut(expression_vector, 
                  breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
                  labels = c("Low", "Medium", "High"),
                  include.lowest = TRUE)
  }
  
  return(as.character(groups))
}

# Function to perform survival analysis for a single chaperone
perform_survival_analysis <- function(project_data, gene_id, clinical_data) {
  project_id <- project_data$project_id
  gene_name <- project_data$gene_names[[gene_id]]
  
  if (is.null(gene_name) || is.na(gene_name)) {
    gene_name <- gene_id
  }
  
  # Get expression data for this gene
  expression_data <- project_data$normalized_expression[gene_id, ]
  
  # Extract patient IDs from expression data (remove sample type suffix)
  sample_ids <- names(expression_data)
  patient_ids <- str_replace(sample_ids, "(-01.+$|-11.+$)", "")
  
  # Create expression dataframe
  expr_df <- data.frame(
    patient_id = patient_ids,
    sample_id = sample_ids,
    expression = as.numeric(expression_data),
    stringsAsFactors = FALSE
  )
  
  # Remove duplicates (keep first occurrence)
  expr_df <- expr_df[!duplicated(expr_df$patient_id), ]
  
  # Clean patient IDs for merging
  expr_df$patient_id_clean <- str_replace(expr_df$patient_id, "^TCGA-", "")
  clinical_data$submitter_id_clean <- str_replace(clinical_data$submitter_id, "^TCGA-", "")
  
  # Merge with clinical data
  merged_data <- merge(expr_df, clinical_data[, c("submitter_id_clean", "vital_status", "days_to_death", "days_to_last_follow_up")], 
                       by.x = "patient_id_clean", by.y = "submitter_id_clean")
  
  # Remove patients with missing expression or survival data
  merged_data <- merged_data[!is.na(merged_data$expression), ]
  
  # Create survival variables
  merged_data$status <- ifelse(merged_data$vital_status == "Dead", 1, 0)
  merged_data$time <- ifelse(!is.na(merged_data$days_to_death), 
                            merged_data$days_to_death, 
                            merged_data$days_to_last_follow_up)
  
  # Remove patients with missing time data
  merged_data <- merged_data[!is.na(merged_data$time) & merged_data$time > 0, ]
  
  if (nrow(merged_data) < 10) {
    message("    Insufficient data for ", gene_name, " in ", project_id)
    return(NULL)
  }
  
  # Create expression groups (using median split)
  merged_data$expression_group <- create_expression_groups(merged_data$expression, gene_name, "median")
  
  # Check if we have both groups
  if (length(unique(merged_data$expression_group)) < 2) {
    message("    Only one expression group for ", gene_name, " in ", project_id)
    return(NULL)
  }
  
  # Perform survival analysis
  tryCatch({
    surv_fit <- survfit(Surv(time, status) ~ expression_group, data = merged_data)
    surv_diff <- survdiff(Surv(time, status) ~ expression_group, data = merged_data)
    p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
    
    # Determine which group has worse survival (lower median survival time)
    # Extract median survival times more robustly
    high_median <- NA
    low_median <- NA
    
    # Try to get median survival from the survival fit object
    tryCatch({
      # Get summary of survival fit
      surv_summary <- summary(surv_fit)
      
      # Extract median survival for each strata
      if (!is.null(surv_fit$table)) {
        median_survival <- surv_fit$table
        
        # Handle matrix case (multiple groups)
        if (is.matrix(median_survival)) {
          if ("median" %in% colnames(median_survival)) {
            strata_names <- rownames(median_survival)
            median_values <- median_survival[, "median"]
            names(median_values) <- strata_names
            
            # Find High and Low groups
            high_idx <- grep("High", strata_names, ignore.case = TRUE)
            low_idx <- grep("Low", strata_names, ignore.case = TRUE)
            
            if (length(high_idx) > 0) {
              high_median <- median_values[high_idx[1]]
            }
            if (length(low_idx) > 0) {
              low_median <- median_values[low_idx[1]]
            }
          }
        } else {
          # Handle vector case (single value)
          if (length(median_survival) == 1 && names(median_survival) == "median") {
            # This shouldn't happen with grouped data, but handle it
            high_median <- median_survival
            low_median <- median_survival
          }
        }
      }
      
      # Alternative approach using quantile function on survival times
      if (is.na(high_median) || is.na(low_median)) {
        # Get the original data and calculate medians manually
        high_patients <- merged_data$time[merged_data$expression_group == "High" & merged_data$status == 1]
        low_patients <- merged_data$time[merged_data$expression_group == "Low" & merged_data$status == 1]
        
        if (length(high_patients) > 0) {
          high_median <- median(high_patients, na.rm = TRUE)
        }
        if (length(low_patients) > 0) {
          low_median <- median(low_patients, na.rm = TRUE)
        }
      }
      
    }, error = function(e) {
      message("      Warning: Could not extract median survival times: ", e$message)
    })
    
    # Handle cases where median is not reached (NA) or not found
    if (is.na(high_median)) high_median <- Inf
    if (is.na(low_median)) low_median <- Inf
    
    # Determine worse survival group
    if (high_median < low_median) {
      worse_survival_group <- "High"
      survival_direction <- "High expression associated with worse survival"
    } else if (low_median < high_median) {
      worse_survival_group <- "Low"
      survival_direction <- "Low expression associated with worse survival"
    } else {
      worse_survival_group <- "Equal"
      survival_direction <- "No clear difference in survival"
    }
    
    # Create survival plot
    p <- ggsurvplot(
      surv_fit,
      data = merged_data,
      pval = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.col = "strata",
      linetype = "strata",
      palette = c("#E41A1C", "#377EB8"),
      title = paste0(project_id, ": ", gene_name, " Expression"),
      xlab = "Days",
      ylab = "Overall Survival Probability",
      legend.title = paste(gene_name, "Expression"),
      legend.labs = c("High", "Low")
    )
    
    # Save plots - PNG only with white background
    png_file <- file.path(individual_dir, paste0(project_id, "_", gene_name, "_survival.png"))
    
    # Save PNG with white background
    png(png_file, width = 10, height = 8, units = "in", res = 300, bg = "white")
    print(p)
    dev.off()
    
    message("    Survival plot created for ", gene_name, " (p = ", round(p_value, 4), ")")
    
    return(list(
      gene_name = gene_name,
      gene_id = gene_id,
      p_value = p_value,
      n_patients = nrow(merged_data),
      worse_survival_group = worse_survival_group,
      survival_direction = survival_direction,
      high_median_survival = if(is.infinite(high_median)) NA else high_median,
      low_median_survival = if(is.infinite(low_median)) NA else low_median,
      plot = p,
      surv_fit = surv_fit
    ))
    
  }, error = function(e) {
    message("    Error in survival analysis for ", gene_name, ": ", e$message)
    return(NULL)
  })
}

# Function to create combined survival plot for all chaperones in a project
create_combined_survival_plot <- function(survival_results, project_id) {
  if (length(survival_results) == 0) {
    return(NULL)
  }
  
  # Create a list of plots
  plot_list <- lapply(survival_results, function(x) {
    if (!is.null(x) && !is.null(x$plot)) {
      # Modify the plot for smaller size in grid
      p <- x$plot$plot + 
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8)
        ) +
        ggtitle(paste(x$gene_name, "\n(p =", round(x$p_value, 4), ")"))
      return(p)
    }
    return(NULL)
  })
  
  # Remove NULL plots
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if (length(plot_list) == 0) {
    return(NULL)
  }
  
  # Calculate grid dimensions
  n_plots <- length(plot_list)
  n_cols <- min(4, ceiling(sqrt(n_plots)))
  n_rows <- ceiling(n_plots / n_cols)
  
  # Create combined plot
  tryCatch({
    # Save combined plot as PNG only with white background
    png(file.path(combined_dir, paste0(project_id, "_all_chaperones_survival.png")),
        width = 4 * n_cols, height = 3 * n_rows + 1, units = "in", res = 300, bg = "white")
    combined_plot <- do.call(grid.arrange, c(plot_list, list(ncol = n_cols)))
    dev.off()
    
    message("  Combined survival plot created for ", project_id, " (", n_plots, " chaperones)")
    return(combined_plot)
    
  }, error = function(e) {
    message("  Error creating combined plot for ", project_id, ": ", e$message)
    return(NULL)
  })
}

# Function to create single overlay plot with all chaperones
create_single_overlay_plot <- function(survival_results, project_id, clinical_data) {
  if (length(survival_results) == 0) {
    return(NULL)
  }
  
  message("  Creating single overlay plot for ", project_id, "...")
  
  # Get project name for title
  project_name <- tryCatch({
    info <- TCGAbiolinks:::getGDCprojects()
    name <- info$name[info$project_id == project_id]
    if(length(name) == 0 || is.na(name)) project_id else name
  }, error = function(e) {
    return(project_id)
  })
  
  # Collect all survival data
  all_surv_data <- data.frame()
  
  # Use a color palette for different chaperones
  n_chaperones <- length(survival_results)
  colors <- rainbow(n_chaperones)
  names(colors) <- names(survival_results)
  
  for (gene_name in names(survival_results)) {
    result <- survival_results[[gene_name]]
    if (is.null(result) || is.null(result$surv_fit)) next
    
    # Extract survival fit data
    surv_fit <- result$surv_fit
    
    # Convert survival fit to data frame
    surv_df <- data.frame(
      time = surv_fit$time,
      n.risk = surv_fit$n.risk,
      n.event = surv_fit$n.event,
      surv = surv_fit$surv,
      strata = rep(names(surv_fit$strata), surv_fit$strata),
      stringsAsFactors = FALSE
    )
    
    # Extract only the "High" expression group
    high_data <- surv_df[grepl("High", surv_df$strata), ]
    
    if (nrow(high_data) > 0) {
      temp_data <- data.frame(
        time = high_data$time,
        survival = high_data$surv,
        gene = gene_name,
        p_value = result$p_value,
        significant = result$p_value < 0.05,
        stringsAsFactors = FALSE
      )
      
      all_surv_data <- rbind(all_surv_data, temp_data)
    }
  }
  
  if (nrow(all_surv_data) == 0) {
    message("  No survival data available for overlay plot")
    return(NULL)
  }
  
  # Create the overlay plot
  tryCatch({
    p <- ggplot(all_surv_data, aes(x = .data$time, y = .data$survival, color = .data$gene)) +
      geom_step(size = 1.2, alpha = 0.8) +
      scale_color_manual(values = colors[unique(all_surv_data$gene)]) +
      labs(
        title = paste0(project_name, ": Chaperone Expression Survival Curves"),
        subtitle = "High expression groups for each chaperone",
        x = "Days",
        y = "Survival Probability",
        color = "Chaperone Gene"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      ylim(0, 1) +
      guides(color = guide_legend(ncol = 1))
    
    # Add significance annotation
    sig_genes <- unique(all_surv_data$gene[all_surv_data$significant])
    if (length(sig_genes) > 0) {
      p <- p + 
        annotate("text", x = Inf, y = 0.05, 
                label = paste("Significant (p<0.05):", paste(sig_genes, collapse = ", ")),
                hjust = 1, vjust = 0, size = 3, color = "red")
    }
    
    # Save overlay plot as PNG only with white background
    ggsave(
      filename = file.path(combined_dir, paste0(project_id, "_chaperones_overlay_survival.png")),
      plot = p,
      width = 12,
      height = 8,
      dpi = 300,
      bg = "white"
    )
    
    message("  Single overlay plot created for ", project_id)
    return(p)
    
  }, error = function(e) {
    message("  Error creating overlay plot for ", project_id, ": ", e$message)
    return(NULL)
  })
}

# Main processing loop
message("\n=== Processing TCGA projects ===")

# Process each TCGA file
all_results <- list()

for (file_path in tcga_files) {
  project_id <- extract_project_id(file_path)
  
  message("\n==== Processing ", project_id, " ====")
  
  # Process the dataset and perform DESeq2 normalization
  project_data <- process_tcga_for_survival(file_path)
  
  if (is.null(project_data)) {
    next
  }
  
  # Get clinical data for the project
  message("  Getting clinical data...")
  clinical_data <- tryCatch({
    GDCquery_clinic(project = project_id, type = "clinical")
  }, error = function(e) {
    message("  Error retrieving clinical data: ", e$message)
    return(NULL)
  })
  
  if (is.null(clinical_data) || nrow(clinical_data) == 0) {
    message("  No clinical data available - skipping")
    next
  }
  
  message("  Found clinical data for ", nrow(clinical_data), " patients")
  
  # Perform survival analysis for each chaperone gene
  survival_results <- list()
  
  for (gene_id in rownames(project_data$normalized_expression)) {
    gene_name <- project_data$gene_names[[gene_id]]
    message("  Analyzing survival for ", gene_name, "...")
    
    result <- perform_survival_analysis(project_data, gene_id, clinical_data)
    if (!is.null(result)) {
      survival_results[[gene_name]] <- result
    }
  }
  
  # Create combined plot for this project
  if (length(survival_results) > 0) {
    create_combined_survival_plot(survival_results, project_id)
    create_single_overlay_plot(survival_results, project_id, clinical_data)
    all_results[[project_id]] <- survival_results
  }
}

# Create summary tables
message("\n=== Creating summary tables ===")

# Table 1: Detailed results for each gene-project combination
detailed_summary <- data.frame()
for (project_id in names(all_results)) {
  for (gene_name in names(all_results[[project_id]])) {
    result <- all_results[[project_id]][[gene_name]]
    detailed_summary <- rbind(detailed_summary, data.frame(
      Project = project_id,
      Gene = gene_name,
      P_Value = result$p_value,
      N_Patients = result$n_patients,
      Significant = result$p_value < 0.05,
      Worse_Survival_Group = result$worse_survival_group,
      Survival_Direction = result$survival_direction,
      High_Median_Survival = result$high_median_survival,
      Low_Median_Survival = result$low_median_survival,
      stringsAsFactors = FALSE
    ))
  }
}

# Table 2: Summary of how many times each chaperone was significant
chaperone_significance_count <- detailed_summary %>%
  group_by(Gene) %>%
  summarise(
    Total_Projects_Tested = n(),
    Times_Significant = sum(Significant),
    Significance_Rate = round(Times_Significant / Total_Projects_Tested * 100, 1),
    Mean_P_Value = round(mean(P_Value, na.rm = TRUE), 4),
    Min_P_Value = round(min(P_Value, na.rm = TRUE), 4),
    Projects_Significant = paste(unique(Project[Significant]), collapse = ", "),
    Times_High_Worse = sum(Significant & Worse_Survival_Group == "High"),
    Times_Low_Worse = sum(Significant & Worse_Survival_Group == "Low"),
    Times_Equal = sum(Significant & Worse_Survival_Group == "Equal"),
    Survival_Direction_Pattern = case_when(
      Times_High_Worse > 0 & Times_Low_Worse == 0 & Times_Equal == 0 ~ "Always High expression worse",
      Times_Low_Worse > 0 & Times_High_Worse == 0 & Times_Equal == 0 ~ "Always Low expression worse", 
      Times_High_Worse > Times_Low_Worse & Times_High_Worse > Times_Equal ~ "Mostly High expression worse",
      Times_Low_Worse > Times_High_Worse & Times_Low_Worse > Times_Equal ~ "Mostly Low expression worse",
      Times_Equal > Times_High_Worse & Times_Equal > Times_Low_Worse ~ "Mostly no clear direction",
      Times_High_Worse == Times_Low_Worse & Times_High_Worse > 0 ~ "Mixed (equal High/Low worse)",
      TRUE ~ "No significant results"
    ),
    Most_Common_Direction = case_when(
      Times_High_Worse > Times_Low_Worse ~ "High expression worse",
      Times_Low_Worse > Times_High_Worse ~ "Low expression worse",
      TRUE ~ "Mixed/Equal"
    ),
    Avg_High_Median_Days = round(mean(High_Median_Survival[Significant & !is.na(High_Median_Survival)], na.rm = TRUE), 0),
    Avg_Low_Median_Days = round(mean(Low_Median_Survival[Significant & !is.na(Low_Median_Survival)], na.rm = TRUE), 0),
    .groups = 'drop'
  ) %>%
  arrange(desc(Times_Significant), Mean_P_Value)

# Save both summary tables
write_csv(detailed_summary, file.path(output_dir, "detailed_survival_analysis_results.csv"))
write_csv(chaperone_significance_count, file.path(output_dir, "chaperone_significance_summary.csv"))

# Create a significance matrix for visualization
significance_matrix <- detailed_summary %>%
  select(Project, Gene, Significant) %>%
  pivot_wider(names_from = Gene, values_from = Significant, values_fill = FALSE) %>%
  column_to_rownames("Project")

# Convert logical to numeric for better visualization
significance_matrix <- as.data.frame(lapply(significance_matrix, as.numeric))

# Save significance matrix
write_csv(cbind(Project = rownames(significance_matrix), significance_matrix), 
          file.path(output_dir, "significance_matrix.csv"))

# Print summary statistics
message("\n=== Summary Statistics ===")
message("Total projects analyzed: ", length(all_results))
message("Total gene-project combinations: ", nrow(detailed_summary))
message("Significant associations (p < 0.05): ", sum(detailed_summary$Significant))
message("Chaperones tested: ", length(unique(detailed_summary$Gene)))

# Print top significant chaperones
message("\n=== Top Significant Chaperones ===")
top_chaperones <- head(chaperone_significance_count, 10)
if (nrow(top_chaperones) > 0) {
  for (i in seq_len(nrow(top_chaperones))) {
    chap <- top_chaperones[i, ]
    message(sprintf("%s: significant in %d/%d projects (%.1f%%) - mean p-value: %.4f", 
                    chap$Gene, chap$Times_Significant, chap$Total_Projects_Tested, 
                    chap$Significance_Rate, chap$Mean_P_Value))
  }
}

# Close logging
sink()
close(con)

message("\nChaperone expression-based survival analysis complete!")
message("Results saved in: ", output_dir)
message("- Individual plots: ", individual_dir)
message("- Combined plots: ", combined_dir)
message("- Detailed results: detailed_survival_analysis_results.csv")
message("- Chaperone significance summary: chaperone_significance_summary.csv")
message("- Significance matrix: significance_matrix.csv")
