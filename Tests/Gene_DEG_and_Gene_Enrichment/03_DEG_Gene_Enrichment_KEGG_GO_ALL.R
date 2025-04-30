library(pathfindR)
library(readr)
library(stringr)

# Base directory for DEG results
base_dir <- "TCGA-Chaperones/Gene_DEG_and_Gene_Enrichment/DEG_results"

# Get all cancer projects (subdirectories)
projects <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
message(paste("Found", length(projects), "projects with DEG results"))

# Track overall progress
total_files <- 0
completed_files <- 0
skipped_files <- 0

# Process each project
for (project in projects) {
  project_path <- file.path(base_dir, project)
  message(paste("Processing project:", project))

  # Get all DEG result files in this project
  deg_files <- list.files(project_path, pattern = "_DEG_results\\.csv$", full.names = TRUE)
  total_files <- total_files + length(deg_files)

  if (length(deg_files) == 0) {
    message(paste("  No DEG result files found in", project))
    next
  }

  # Process each gene's DEG file
  for (file_path in deg_files) {
    # Extract gene name from file name
    file_name <- basename(file_path)
    gene_name <- gsub("_DEG_results\\.csv$", "", file_name)

    # Check if ANY pathfindR result folder exists for this gene (including numbered ones)
    pathfindR_folder_pattern <- paste0("^", gene_name, "_pathfindR(\\([0-9]+\\))?$")
    existing_folders <- list.dirs(project_path, full.names = FALSE, recursive = FALSE)
    matching_folders <- grep(pathfindR_folder_pattern, existing_folders, value = TRUE)

    # Check if ANY of those folders contain a results file
    has_results <- FALSE
    for (folder in matching_folders) {
      potential_result_file <- file.path(project_path, folder, paste0(gene_name, "_pathfindR_results.csv"))
      if (file.exists(potential_result_file)) {
        has_results <- TRUE
        break
      }
    }

    # Skip if we already have results
    if (has_results) {
      message(paste("  Skipping gene:", gene_name, "(results already exist)"))
      skipped_files <- skipped_files + 1
      next
    }

    # If we have partial folders without results, we should remove them
    # This is to prevent pathfindR from creating too many duplicate folders
    if (length(matching_folders) > 0) {
      message(paste("  Cleaning up incomplete folders for:", gene_name))
      for (folder in matching_folders) {
        unlink(file.path(project_path, folder), recursive = TRUE)
      }
    }

    # Create a subdirectory for this gene's pathfindR results
    pathfindR_dir <- file.path(project_path, paste0(gene_name, "_pathfindR"))
    result_file <- file.path(pathfindR_dir, paste0(gene_name, "_pathfindR_results.csv"))

    message(paste("  Processing gene:", gene_name))

    # Read the DEG results
    tryCatch(
      {
        deg_data <- read_csv(file_path, show_col_types = FALSE)

        # Check if we have the required columns
        if (all(c("gene_name", "logFC", "FDR") %in% colnames(deg_data))) {
          # Create input for pathfindR
          input_df <- data.frame(
            Gene_symbol = deg_data$gene_name,
            logFC = deg_data$logFC,
            FDR_p = deg_data$FDR
          )
        } else if (all(c("gene_name", "logFC", "PValue") %in% colnames(deg_data))) {
          # Create input for pathfindR
          input_df <- data.frame(
            Gene_symbol = deg_data$gene_name,
            logFC = deg_data$logFC,
            FDR_p = deg_data$PValue
          )
        } else {
          message("    Required columns not found in file, skipping...")
          next
        }

        # Remove rows with NA in Gene_symbol
        input_df <- input_df[!is.na(input_df$Gene_symbol), ]

        # Check if there are enough genes for analysis
        if (nrow(input_df) < 10) {
          message("    Too few genes with valid symbols (", nrow(input_df), "), skipping...")
          next
        }

        # Ensure the directory is fresh for this run
        if (dir.exists(pathfindR_dir)) {
          unlink(pathfindR_dir, recursive = TRUE)
        }
        dir.create(pathfindR_dir, showWarnings = FALSE)

        message("    Running pathfindR analysis with GO-All and KEGG PIN...")

        # Run pathfindR analysis
        tryCatch(
          {
            # Run the analysis with specified parameters
            pathfinder_results <- run_pathfindR(
              input_df,
              output_dir = pathfindR_dir,
              gene_sets = "GO-All",
              pin_name = "KEGG",
              min_gset_size = 10,
              max_gset_size = 300,
              p_val_threshold = 0.05,
              enrichment_threshold = 0, # Disables term-gene heatmaps
              silent_option = TRUE # Suppresses other visualizations
            )

            # Find the actual folder where results were saved
            # (which might be pathfindR_dir or pathfindR_dir(1), etc.)
            result_folders <- list.dirs(project_path, full.names = FALSE, recursive = FALSE)
            result_folder_pattern <- paste0("^", gene_name, "_pathfindR(\\([0-9]+\\))?$")
            matching_result_folders <- grep(result_folder_pattern, result_folders, value = TRUE)

            if (length(matching_result_folders) > 0) {
              # Get the most recently created folder
              folder_stats <- file.info(file.path(project_path, matching_result_folders))
              newest_folder <- rownames(folder_stats)[which.max(folder_stats$mtime)]
              newest_folder_name <- basename(newest_folder)

              # Check if we have results to save
              if (nrow(pathfinder_results) > 0) {
                # Create the results filename in the correct folder
                final_result_file <- file.path(
                  project_path, newest_folder_name,
                  paste0(gene_name, "_pathfindR_results.csv")
                )
                write_csv(pathfinder_results, final_result_file)
                message(paste(
                  "    Saved", nrow(pathfinder_results), "enriched terms to",
                  basename(final_result_file), "in folder", newest_folder_name
                ))
                completed_files <- completed_files + 1
              } else {
                message("    No enriched terms found")
                # Create an empty file as a marker that we processed this gene
                final_result_file <- file.path(
                  project_path, newest_folder_name,
                  paste0(gene_name, "_pathfindR_results.csv")
                )
                write_csv(data.frame(), final_result_file)
                completed_files <- completed_files + 1
              }
            } else {
              message("    Warning: Could not find results folder!")
            }
          },
          error = function(e) {
            message(paste("    Error in pathfindR analysis:", e$message))
          }
        )
      },
      error = function(e) {
        message(paste("    Error reading file:", e$message))
      }
    )
  }

  message(paste("Completed analysis for project:", project))
}

message(paste("All pathfindR enrichment analyses complete!"))
message(paste("Processed:", completed_files, "files"))
message(paste("Skipped (already existed):", skipped_files, "files"))
message(paste("Total files:", total_files))
