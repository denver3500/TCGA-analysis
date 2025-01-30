library(SummarizedExperiment)
library(TCGAbiolinks)

# Get the list of project IDs
Links <- TCGAbiolinks:::getGDCprojects()$project_id

# Create an empty list to store queries
queries <- list()

# Loop through each project in Links to query and download data
for (project in Links) {
  message(paste("Processing project:", project))
  
  # Create a query for the current project
  query <- tryCatch({
    GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal")
    )
  }, error = function(e) {
    message(paste("Error creating query for project:", project))
    message(e)
    return(NULL)
  })
  
  if (!is.null(query)) {
    # Store the query in the list
    queries[[project]] <- query
    
    # Download the data for the current project
    tryCatch({
      GDCdownload(
        query = query,
        files.per.chunk = 100
      )
    }, error = function(e) {
      message(paste("Error downloading data for project:", project))
      message(e)query
    })
  }
}

# Loop through each project in Links to prepare and save data
for (project in Links) {
  if (!is.null(queries[[project]])) {
    # Prepare and save the data for the current project
    tryCatch({
      transcriptomic.exp <- GDCprepare(
        query = queries[[project]]
      )
      # Save the data (example: save as RDS file)
      saveRDS(transcriptomic.exp, file = paste0(project, "_transcriptomic_exp.rds"))
    }, error = function(e) {
      message(paste("Error preparing data for project:", project))
      message(e)
    })
  }
}