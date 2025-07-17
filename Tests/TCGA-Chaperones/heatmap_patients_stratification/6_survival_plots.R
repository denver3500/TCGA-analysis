# Load required libraries
library(TCGAbiolinks)
library(survival)
library(dplyr)
library(readr)
library(stringr)

# Set working directory
setwd("/media/windows/BIOINFORMATICS/TCGA/TCGA-analysis/TCGA-Chaperones/heatmap_patients_stratification")

# Create output directory
output_dir <- "survival_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get list of projects with clustering data
cluster_dir <- "deseq2_clustering"
projects <- list.dirs(cluster_dir, full.names = FALSE, recursive = FALSE)
projects <- projects[projects != ""]

cat("Found", length(projects), "TCGA projects with clustering results\n")

# Process each project
for(project in projects) {
  cat("\n==== Processing", project, "====\n")
  
  # 1. GET CLUSTER DATA
  # Check if cluster file exists
  cluster_file <- file.path(cluster_dir, project, paste0(project, "_patient_clusters.csv"))
  
  if(!file.exists(cluster_file)) {
    cat("  No cluster file found for", project, "\n")
    next
  }
  
  # Read clustering data
  cluster_data <- read_csv(cluster_file, show_col_types = FALSE)
  cat("  Found", nrow(cluster_data), "patients in clusters\n")
  
  # 2. TRIM NAMES AND REMOVE DUPLICATES
  # Clean up patient IDs (remove sample-specific suffix)
  cluster_data$Patient_ID <- str_replace(cluster_data$Patient_ID, "(-01.+$|-11.+$)", "")
  
  # Remove duplicate patient entries
  if(any(duplicated(cluster_data$Patient_ID))) {
    cat("  Removing", sum(duplicated(cluster_data$Patient_ID)), "duplicate patients\n")
    cluster_data <- cluster_data[!duplicated(cluster_data$Patient_ID),]
  }
  
  # Get clinical data for the project
  cat("  Getting clinical data...\n")
  clin <- tryCatch({
    GDCquery_clinic(project = project, type = "clinical")
  }, error = function(e) {
    cat("  Error retrieving clinical data:", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(clin) || nrow(clin) == 0) {
    cat("  No clinical data available - skipping\n")
    next
  }
  
  # Prepare patient groups for survival analysis
  # Clean up patient IDs for joining
  cluster_data$Patient_ID_clean <- str_replace(cluster_data$Patient_ID, "^TCGA-", "")
  clin$submitter_id_clean <- str_replace(clin$submitter_id, "^TCGA-", "")
  
  # Create groups list for TCGAanalyze_survival
  groups <- list()
  for(cluster_num in sort(unique(cluster_data$Patient_Cluster))) {
    # Get TCGA patient barcodes for this cluster
    cluster_patients <- cluster_data$Patient_ID[cluster_data$Patient_Cluster == cluster_num]
    
    # Make sure all barcodes start with "TCGA-"
    cluster_patients <- sapply(cluster_patients, function(id) {
      if(!startsWith(id, "TCGA-")) return(paste0("TCGA-", id))
      return(id)
    })
    
    # Add to groups list
    groups[[paste0("Cluster_", cluster_num)]] <- unique(cluster_patients)
  }
  
  # Check group sizes but don't skip - just warn
  group_sizes <- sapply(groups, length)
  min_group_size <- min(group_sizes)
  if(min_group_size < 5) {
    cat("  Warning: Some clusters have fewer than 5 patients (min:", min_group_size, ")\n")
    cat("  Group sizes:", paste(names(group_sizes), "=", group_sizes, collapse=", "), "\n")
  }
  
  # Skip only if any group is completely empty
  if(min_group_size == 0) {
    cat("  Skipping - found empty clusters\n")
    next
  }
  
  # 3. RUN SURVIVAL ANALYSIS
  # Create an empty data frame for TCGAanalyze_survival's data parameter
  cat("  Creating survival plot...\n")
  all_patients <- unique(unlist(groups))
  dummy_data <- data.frame(row.names = all_patients)
  
  # Colors for the clusters
  cluster_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                      "#FFFF33", "#A65628", "#F781BF")[1:length(groups)]
  names(cluster_colors) <- names(groups)
  
  # Get project name
  project_name <- tryCatch({
    info <- TCGAbiolinks:::getGDCprojects()
    name <- info$name[info$project_id == project]
    if(length(name) == 0 || is.na(name)) project else name
  }, error = function(e) {
    return(project)
  })
  
  # Run TCGAanalyze_survival for overall survival
  tryCatch({
    TCGAanalyze_survival(
      clinical_patient = clin,
      groups = groups,
      data = dummy_data,
      main = paste0(project_name, ": Overall Survival by Patient Cluster"),
      filename = file.path(output_dir, paste0(project, "_OS_survival.pdf")),
      color = cluster_colors,
      height = 10,
      width = 12,
      legend = paste("Cluster", sort(unique(cluster_data$Patient_Cluster))),
      is_OS = TRUE
    )
    
    # 4. SAVE IMAGE (also as PNG)
    cat("  Saving PNG version...\n")
    png_file <- file.path(output_dir, paste0(project, "_OS_survival.png"))
    
    # Re-run to generate PNG
    TCGAanalyze_survival(
      clinical_patient = clin,
      groups = groups,
      data = dummy_data,
      main = paste0(project_name, ": Overall Survival by Patient Cluster"),
      filename = png_file,
      color = cluster_colors,
      height = 10,
      width = 12,
      legend = paste("Cluster", sort(unique(cluster_data$Patient_Cluster))),
      is_OS = TRUE,
      imgtype = "png",
      dpi = 300
    )
    
    cat("  Success! Survival plot saved for", project, "\n")
    
    # If the project has progression-free interval data, create that plot too
    if("PFI" %in% colnames(clin) && "PFI.time" %in% colnames(clin)) {
      TCGAanalyze_survival(
        clinical_patient = clin,
        groups = groups,
        data = dummy_data,
        main = paste0(project_name, ": Progression-Free Survival by Patient Cluster"),
        filename = file.path(output_dir, paste0(project, "_PFI_survival.pdf")),
        color = cluster_colors,
        height = 10,
        width = 12,
        legend = paste("Cluster", sort(unique(cluster_data$Patient_Cluster))),
        is_OS = FALSE
      )
      
      # PNG version
      TCGAanalyze_survival(
        clinical_patient = clin,
        groups = groups,
        data = dummy_data,
        main = paste0(project_name, ": Progression-Free Survival by Patient Cluster"),
        filename = file.path(output_dir, paste0(project, "_PFI_survival.png")),
        color = cluster_colors,
        height = 10,
        width = 12,
        legend = paste("Cluster", sort(unique(cluster_data$Patient_Cluster))),
        is_OS = FALSE,
        imgtype = "png",
        dpi = 300
      )
      cat("  Also saved progression-free survival plot\n")
    }
    
  }, error = function(e) {
    cat("  Error in survival analysis:", e$message, "\n")
    
    # Try directly using survival package if TCGAanalyze_survival fails
    cat("  Attempting to create survival plot using base survival functions...\n")
    
    # Clean and merge data
    clin_subset <- clin[, c("submitter_id_clean", "vital_status", "days_to_death", "days_to_last_follow_up")]
    merged <- merge(cluster_data, clin_subset, by.x = "Patient_ID_clean", by.y = "submitter_id_clean")
    
    # Create survival object
    merged$status <- ifelse(merged$vital_status == "Dead", 1, 0)
    merged$time <- ifelse(!is.na(merged$days_to_death), merged$days_to_death, merged$days_to_last_follow_up)
    
    # Skip if missing data
    if(sum(is.na(merged$time) | is.na(merged$status)) > nrow(merged)/2) {
      cat("  Too much missing data for survival analysis\n")
      return()
    }
    
    # Create survival plot
    pdf(file.path(output_dir, paste0(project, "_OS_survival_basic.pdf")), width=10, height=8)
    surv_fit <- survfit(Surv(time, status) ~ Patient_Cluster, data=merged)
    plot(surv_fit, col=1:length(unique(merged$Patient_Cluster)), lwd=2,
         main=paste0(project_name, ": Overall Survival by Patient Cluster"),
         xlab="Days", ylab="Survival Probability")
    legend("topright", legend=paste("Cluster", sort(unique(merged$Patient_Cluster))), 
           col=1:length(unique(merged$Patient_Cluster)), lwd=2)
    dev.off()
    
    cat("  Created basic survival plot instead\n")
  })
}

cat("\nSurvival analysis complete! Results are saved in the", output_dir, "directory.\n")