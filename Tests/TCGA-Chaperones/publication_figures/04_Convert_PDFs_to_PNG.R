#!/usr/bin/env Rscript

# Convert PDF heatmaps to PNG for easier integration into publication figures

library(magick)

# Define paths
stratification_base <- "TCGA-Chaperones/heatmap_patients_stratification/deseq2_clustering"
pathway_base <- "TCGA-Chaperones/heatmap_patients_stratification"

message("Converting PDF heatmaps to PNG for publication figures...")

# Convert patient stratification heatmaps
cancer_types <- list.dirs(stratification_base, full.names = FALSE, recursive = FALSE)
cancer_types <- cancer_types[grepl("^TCGA-", cancer_types)]

for (cancer in cancer_types) {
  pdf_file <- file.path(stratification_base, cancer, 
                       paste0(cancer, "_bidirectional_clustered_heatmap.pdf"))
  png_file <- file.path(stratification_base, cancer, 
                       paste0(cancer, "_bidirectional_clustered_heatmap.png"))
  
  if (file.exists(pdf_file) && !file.exists(png_file)) {
    tryCatch({
      img <- image_read_pdf(pdf_file, density = 300)
      image_write(img, png_file, format = "png")
      message("Converted: ", cancer, " stratification heatmap")
    }, error = function(e) {
      message("Error converting ", cancer, ": ", e$message)
    })
  }
}

# Convert pathway heatmaps
pathway_files <- c("pathways_heatmap.pdf", "top_enriched_pathways.pdf")

for (pdf_name in pathway_files) {
  pdf_file <- file.path(pathway_base, pdf_name)
  png_file <- file.path(pathway_base, gsub("\\.pdf$", ".png", pdf_name))
  
  if (file.exists(pdf_file) && !file.exists(png_file)) {
    tryCatch({
      img <- image_read_pdf(pdf_file, density = 300)
      image_write(img, png_file, format = "png")
      message("Converted: ", pdf_name)
    }, error = function(e) {
      message("Error converting ", pdf_name, ": ", e$message)
    })
  }
}

message("PDF to PNG conversion completed!")
