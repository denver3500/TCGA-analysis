library(SummarizedExperiment)
library(TCGAbiolinks)

# Get the list of project IDs
Links <- TCGAbiolinks:::getGDCprojects()$project_id

# Create an empty list to store queries
queries <- list()

# Loop through each project in Links to query and download data
for (project in Links) {
  # Create a query for the current project
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  
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
    message(e)
  })
}

# Loop through each project in Links to prepare and save data
for (project in Links) {
  # Prepare and save the data for the current project
  tryCatch({
    transcriptomic.exp <- GDCprepare(
      query = queries[[project]],
      save = TRUE,
      save.filename = paste0(project, "_transcriptomicexp.rda")
    )
  }, error = function(e) {
    message(paste("Error preparing data for project:", project))
    message(e)
  })
}