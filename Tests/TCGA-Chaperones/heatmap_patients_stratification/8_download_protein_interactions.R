library(biomaRt)
library(readr)
library(dplyr)

# Read the TSV file
interactions <- read_tsv("TCGA-Chaperones/heatmap_patients_stratification/chaperone_interactors_raw.tsv", 
                        show_col_types = FALSE)

# Step 1: Select only the two columns we need
interactions_clean <- interactions %>%
  select(`# ID(s) interactor A`, `ID(s) interactor B`) %>%
  rename(interactor_A = `# ID(s) interactor A`,
         interactor_B = `ID(s) interactor B`)

# Step 2: Filter rows that contain "uniprotkb:" in both columns
interactions_filtered <- interactions_clean %>%
  filter(grepl("uniprotkb:", interactor_A) & grepl("uniprotkb:", interactor_B))

# Step 3: Extract UniProt IDs (remove "uniprotkb:" prefix)
extract_uniprot_id <- function(x) {
  # Extract the part after "uniprotkb:"
  uniprot_part <- sub(".*uniprotkb:([^|]*)", "\\1", x)
  return(uniprot_part)
}

interactions_trimmed <- interactions_filtered %>%
  mutate(
    uniprot_A = extract_uniprot_id(interactor_A),
    uniprot_B = extract_uniprot_id(interactor_B)
  ) %>%
  select(uniprot_A, uniprot_B)

# Step 4: Convert UniProt IDs to Ensembl Gene IDs using biomaRt
message("Connecting to biomaRt...")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get all unique UniProt IDs
all_uniprot_ids <- unique(c(interactions_trimmed$uniprot_A, interactions_trimmed$uniprot_B))
message("Converting ", length(all_uniprot_ids), " unique UniProt IDs to Ensembl...")

# Query biomaRt to get Ensembl gene IDs
uniprot_to_ensembl <- getBM(
  attributes = c("uniprotswissprot", "ensembl_gene_id"),
  filters = "uniprotswissprot",
  values = all_uniprot_ids,
  mart = mart
)

# Create a mapping function
map_to_ensembl <- function(uniprot_id) {
  ensembl_id <- uniprot_to_ensembl$ensembl_gene_id[uniprot_to_ensembl$uniprotswissprot == uniprot_id]
  if(length(ensembl_id) > 0) {
    return(ensembl_id[1])  # Take first match if multiple
  } else {
    return(NA)
  }
}

# Add Ensembl columns
interactions_with_ensembl <- interactions_trimmed %>%
  mutate(
    ensembl_A = sapply(uniprot_A, map_to_ensembl),
    ensembl_B = sapply(uniprot_B, map_to_ensembl)
  ) %>%
  # Remove rows where we couldn't map either protein
  filter(!is.na(ensembl_A) & !is.na(ensembl_B))

# Step 5: Load your chaperone gene list to identify which are chaperones
chaperone_genes <- read_csv("TCGA-Chaperones/gene_list.csv", show_col_types = FALSE)
chaperone_ensembl_ids <- chaperone_genes$ENSEMBL

# Create final dataframe with chaperone vs interactor organization
organize_interactions <- function(ensembl_a, ensembl_b, uniprot_a, uniprot_b) {
  # Check which one is a chaperone
  a_is_chaperone <- ensembl_a %in% chaperone_ensembl_ids
  b_is_chaperone <- ensembl_b %in% chaperone_ensembl_ids
  
  if(a_is_chaperone && !b_is_chaperone) {
    return(data.frame(
      chaperone_ensembl = ensembl_a,
      chaperone_uniprot = uniprot_a,
      interactor_ensembl = ensembl_b,
      interactor_uniprot = uniprot_b
    ))
  } else if(!a_is_chaperone && b_is_chaperone) {
    return(data.frame(
      chaperone_ensembl = ensembl_b,
      chaperone_uniprot = uniprot_b,
      interactor_ensembl = ensembl_a,
      interactor_uniprot = uniprot_a
    ))
  } else if(a_is_chaperone && b_is_chaperone) {
    # Both are chaperones - keep both directions
    return(data.frame(
      chaperone_ensembl = c(ensembl_a, ensembl_b),
      chaperone_uniprot = c(uniprot_a, uniprot_b),
      interactor_ensembl = c(ensembl_b, ensembl_a),
      interactor_uniprot = c(uniprot_b, uniprot_a)
    ))
  } else {
    # Neither is a chaperone - skip this interaction
    return(NULL)
  }
}

# Apply the organization function
final_interactions <- NULL
for(i in 1:nrow(interactions_with_ensembl)) {
  row_result <- organize_interactions(
    interactions_with_ensembl$ensembl_A[i],
    interactions_with_ensembl$ensembl_B[i],
    interactions_with_ensembl$uniprot_A[i],
    interactions_with_ensembl$uniprot_B[i]
  )
  if(!is.null(row_result)) {
    final_interactions <- rbind(final_interactions, row_result)
  }
}

# Remove duplicates
final_interactions <- final_interactions %>%
  distinct(chaperone_ensembl, interactor_ensembl, .keep_all = TRUE)

# Add gene names if available
message("Adding gene names...")
all_ensembl_for_names <- unique(c(final_interactions$chaperone_ensembl, 
                                 final_interactions$interactor_ensembl))

ensembl_to_names <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = all_ensembl_for_names,
  mart = mart
)

# Add gene names
final_interactions <- final_interactions %>%
  left_join(ensembl_to_names, by = c("chaperone_ensembl" = "ensembl_gene_id")) %>%
  rename(chaperone_gene_name = external_gene_name) %>%
  left_join(ensembl_to_names, by = c("interactor_ensembl" = "ensembl_gene_id")) %>%
  rename(interactor_gene_name = external_gene_name)

# Reorder columns for better readability
final_interactions <- final_interactions %>%
  select(chaperone_ensembl, chaperone_gene_name, chaperone_uniprot,
         interactor_ensembl, interactor_gene_name, interactor_uniprot)

# Save the results
write_csv(final_interactions, "TCGA-Chaperones/heatmap_patients_stratification/chaperone_protein_interactions.csv")

# Print summary
message("Processing complete!")
message("Original interactions: ", nrow(interactions_clean))
message("After filtering for UniProt IDs: ", nrow(interactions_filtered))
message("After Ensembl mapping: ", nrow(interactions_with_ensembl))
message("Final chaperone interactions: ", nrow(final_interactions))
message("Unique chaperones involved: ", length(unique(final_interactions$chaperone_ensembl)))
message("Unique interactors: ", length(unique(final_interactions$interactor_ensembl)))

# Show preview
message("\nPreview of results:")
print(head(final_interactions, 10))

# Save summary statistics
summary_stats <- data.frame(
  metric = c("Original_interactions", "UniProt_filtered", "Ensembl_mapped", 
            "Final_chaperone_interactions", "Unique_chaperones", "Unique_interactors"),
  count = c(nrow(interactions_clean), nrow(interactions_filtered), 
           nrow(interactions_with_ensembl), nrow(final_interactions),
           length(unique(final_interactions$chaperone_ensembl)),
           length(unique(final_interactions$interactor_ensembl)))
)

write_csv(summary_stats, "TCGA-Chaperones/heatmap_patients_stratification/interaction_processing_summary.csv")

message("\nFiles saved:")
message("- TCGA-Chaperones/chaperone_protein_interactions.csv")
message("- TCGA-Chaperones/interaction_processing_summary.csv")