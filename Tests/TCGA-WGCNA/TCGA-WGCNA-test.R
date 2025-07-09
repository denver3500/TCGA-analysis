library(tidyverse)
library(magrittr)
library(WGCNA)
library(DESeq2)
library(genefilter)


data <- readRDS(file.path("TCGA-WGCNA/split_datasets","TCGA-KICH_transcriptomic_exp_processed.rds"))
assayNames(data)[assayNames(data) == "unstranded"] <- "counts"

# Check sample types available (optional, to verify)
table(data$sample_type)

# Filter to keep only Primary Tumor samples
data <- data[, data$sample_type == "Primary Tumor"]

# Create DESeq dataset with filtered data (design=~1 since we now have only one group)
dds <- DESeqDataSet(data, design = ~ 1)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)

# Continue with the rest of your analysis as before
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75) 
q95_wpn <- quantile( rowVars(wpn_vsd), .95)
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

input_mat = t(expr_normalized)
allowWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]
#>            gene_id    colors
#> 1 AC149818.2_FG001      blue
#> 2 AC149829.2_FG003      blue
#> 3 AC182617.3_FG001      blue
#> 4 AC186512.3_FG001 turquoise
#> 5 AC186512.3_FG007 turquoise

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")