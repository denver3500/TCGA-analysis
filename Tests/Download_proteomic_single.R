################################# Libraries
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)

################################# Download Proteome data

query.proteomic.chol <- GDCquery(
  project = "TCGA-CHOL", #Adenomas and Adenocarcinomas
  data.category = "Proteome Profiling", #Proteome data
  data.type = "Protein Expression Quantification",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)

GDCdownload(
  query = query.proteomic.chol,
  files.per.chunk = 100
)

############################### Analysis of Proteome data

chol.proteomic.exp <- GDCprepare(
  query = query.proteomic.chol, 
  save = TRUE, 
  save.filename = "choltranscriptomicexp.rda"
)

protein.id <- read.csv("ProteinID.csv")
id_column <- "AGID"
filtered_data <- chol.proteomic.exp %>%
  filter(AGID %in% protein.id[[id_column]])

write.csv(filtered_data, "Filtered_Proteome_Data.csv", 
          row.names = FALSE)
#########################################################