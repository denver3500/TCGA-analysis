library(SummarizedExperiment)
library(TCGAbiolinks)

#https://docs.gdc.cancer.gov/Data_Portal/PDF/Data_Portal_UG.pdf docs for GDC.
#https://bioconductor.org/packages/release/bioc/manuals/TCGAbiolinks/man/TCGAbiolinks.pdf docs for TCGAbiolinks package
#https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/index.html vignettes for TCGAbiolinks package

query.brca <- GDCquery(
  project = "TCGA-BRCA", #Breast Invasive Carcinoma
  data.category = "Transcriptome Profiling", #Data filter (Methylation/Proteome/Transcriptome/etc)
  data.type = "Gene Expression Quantification", #Data filter (Slide images/Somatic mutations/etc)
  workflow.type = "STAR - Counts",#STAR - Counts is a pipeline used for to quantify gene expression from RNA-Seq data. Quantification is performed with STAR using tumor or normal alignments and generates transcriptome profiling data.
  sample.type = c("Primary Tumor","Solid Tissue Normal") #Sample Type has been removed from the Associated Cases/Biospecimens table and replaced with Tissue Type and Tumor descriptor???
)
GDCdownload(
  query = query.brca,
  files.per.chunk = 100
)

brca.exp <- GDCprepare(
  query = query.brca, 
  save = TRUE, 
  save.filename = "brcaExp.rda"
)

# get subtype information
infomation.subtype <- TCGAquery_subtype(tumor = "BRCA")

# get clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-BRCA",type = "clinical") 

# Which samples are Primary Tumor
samples.primary.tumour <- brca.exp$barcode[brca.exp$shortLetterCode == "TP"]

# which samples are solid tissue normal
samples.solid.tissue.normal <- brca.exp$barcode[brca.exp$shortLetterCode == "NT"]