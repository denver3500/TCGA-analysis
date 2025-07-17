# WGCNA by Patient Clusters - Test Analysis (TCGA-BRCA & TCGA-BLCA)

This directory contains scripts for testing the WGCNA cluster-based analysis pipeline on two specific cancer types: **TCGA-BRCA** (Breast Cancer) and **TCGA-BLCA** (Bladder Cancer).

## Test Directory Structure
```
WGCNA_by_clusters_test/
├── 01_Normalization_by_clusters_test.R       # Data normalization for BRCA & BLCA
├── 02_Network_construction_test.R            # WGCNA network construction for test projects
├── 03_Chaperone_module_analysis_test.R       # Chaperone analysis for test projects
├── processed_data/                           # Normalized data for each cluster
├── results/                                  # WGCNA results and networks
└── chaperone_analysis/                       # Chaperone-specific analysis
```

## Quick Test Pipeline

This is a streamlined version of the main pipeline designed for rapid testing and validation with two representative cancer types.

### Test Projects:
- **TCGA-BRCA**: Breast Invasive Carcinoma
- **TCGA-BLCA**: Bladder Urothelial Carcinoma

### Expected Outputs:

#### For TCGA-BRCA:
- **Cluster 1 Network**: Co-expression modules for BRCA patients in cluster 1
- **Cluster 2 Network**: Co-expression modules for BRCA patients in cluster 2
- **Chaperone Analysis**: How chaperone genes cluster in each patient subgroup

#### For TCGA-BLCA:
- **Cluster 1 Network**: Co-expression modules for BLCA patients in cluster 1  
- **Cluster 2 Network**: Co-expression modules for BLCA patients in cluster 2
- **Chaperone Analysis**: How chaperone genes cluster in each patient subgroup

### Quick Usage:

```bash
# Navigate to test directory
cd TCGA-Chaperones/WGCNA_by_clusters_test

# Run test pipeline (in order):
Rscript 01_Normalization_by_clusters_test.R
Rscript 02_Network_construction_test.R
Rscript 03_Chaperone_module_analysis_test.R
```

### Key Features:

1. **Limited Scope**: Only processes TCGA-BRCA and TCGA-BLCA for faster execution
2. **Same Methodology**: Uses identical analysis methods as the full pipeline
3. **Quick Validation**: Allows you to verify the approach works before running on all cancer types
4. **Representative Results**: Breast and bladder cancers provide good biological diversity

### Expected Runtime:
- **Step 1 (Normalization)**: ~5-10 minutes per project
- **Step 2 (Network Construction)**: ~10-20 minutes per cluster
- **Step 3 (Chaperone Analysis)**: ~2-5 minutes per cluster

### Test Results to Expect:

#### Example for TCGA-BLCA (from your cluster file):
- **Cluster 1**: ~198 patients
- **Cluster 2**: ~209 patients  
- **Networks**: Separate co-expression networks for each cluster
- **Chaperone Modules**: Distribution of your 19 chaperone genes across modules

#### Biological Questions Answered:
1. **Do patient clusters have different co-expression patterns?**
2. **Are chaperone genes consistently clustered together?**
3. **Do chaperone clustering patterns differ between patient subgroups?**
4. **Are patterns similar between breast and bladder cancers?**

### Success Criteria:
✅ **Both projects process successfully**  
✅ **Networks constructed for each cluster**  
✅ **Chaperone genes found in modules**  
✅ **Clear visualizations generated**  
✅ **Summary statistics look reasonable**  

### Next Steps:
If test results look good, you can then run the full pipeline on all cancer types using the main `WGCNA_by_clusters/` directory.

### Troubleshooting:
- **Missing cluster files**: Check if patient cluster CSV files exist
- **No RDS data**: Verify RDS files are available for both projects
- **Memory issues**: Reduce the number of genes or samples if needed
- **WGCNA errors**: Check soft-thresholding power selection

This test setup allows you to validate the entire analytical approach quickly before committing to the full analysis across all TCGA cancer types.
