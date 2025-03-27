################################# Libraries
library(SummarizedExperiment)
library(TCGAbiolinks)

################################# All projects IDs
projects_info <- TCGAbiolinks:::getGDCprojects()
project_ids <- projects_info$project_id
project_ids <- projects_info$project_id[startsWith(projects_info$project_id, "TCGA")] #Filter to leave TCGA data only

#################################  Loop through each project to query and download data
for (project_id in project_ids) {
  message(paste("Processing project:", project_id)) # Print the current project ID
  
  query <- tryCatch({ #
    GDCquery(
      project = project_id,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal") # Include only primary tumor and normal samples
    )
  }, error = function(e) {
    message(paste("Error creating query for project:", project_id)) # Print error message if there's no valid data
    message(e)
    return(NULL)
  })
  
  if (!is.null(query)) { # Check if the query was successful
   
    tryCatch({
      GDCdownload(
        directory = "raw_data/GDCdata",
        query = query,
        files.per.chunk = 100
      )
    }, error = function(e) {
      message(paste("Error downloading data for project:", project_id))
      message(e)
    })
    
    tryCatch({
      transcriptomic_exp <- 
      GDCprepare(
        directory = "raw_data/GDCdata",
        query = query)
      
      # Save the data as RDS file
      saveRDS(transcriptomic_exp, file = file.path("TCGA-Chaperones/rds", paste0(project_id, "_transcriptomic_exp.rds")))
   }, error = function(e) {
      message(paste("Error preparing data for project:", project_id))
      message(e)
    })
  }

}