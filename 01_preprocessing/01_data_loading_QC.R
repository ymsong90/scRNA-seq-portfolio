################################################################################
# Single-Cell RNA-seq Analysis Pipeline
# Step 01: Data Loading and Quality Control
# 
# Description:
#   - Load 10X Genomics data for WT and KO conditions
#   - Perform quality control filtering
#   - Generate QC metrics and visualizations
#
# Input:
#   - ./data/wt/filtered_feature_bc_matrix/
#   - ./data/ko/filtered_feature_bc_matrix/
#
# Output:
#   - porcn.combined: Merged Seurat object (before QC)
#   - QC plots: Violin plots and scatter plots
#
# Author: YMS
# Date: 2025
################################################################################

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(hdf5r)
})

# Set parameters ---------------------------------------------------------------
QC_PARAMS <- list(
    min_features = 200,      # Minimum genes per cell
    max_features = 8000,     # Maximum genes per cell
    max_mt_pct   = 20        # Maximum mitochondrial percentage
)

# Create output directory
if (!dir.exists("./results/01_QC")) dir.create("./results/01_QC", recursive = TRUE)

################################################################################
# 1. Load 10X Data
################################################################################

cat("\n=== Step 1: Loading 10X Data ===\n")

# Load WT data
porcn_wt <- Read10X(data.dir = "./data/wt/filtered_feature_bc_matrix/")
cat("✓ WT data loaded:", ncol(porcn_wt), "cells,", nrow(porcn_wt), "genes\n")

# Load KO data
porcn_ko <- Read10X(data.dir = "./data/ko/filtered_feature_bc_matrix/")
cat("✓ KO data loaded:", ncol(porcn_ko), "cells,", nrow(porcn_ko), "genes\n")

################################################################################
# 2. Create Seurat Objects
################################################################################

cat("\n=== Step 2: Creating Seurat Objects ===\n")

# Create Seurat objects with initial filtering
porcn_wt <- CreateSeuratObject(
    counts      = porcn_wt, 
    project     = 'porcn_wt',
    min.cells   = 3,         # Gene must be detected in at least 3 cells
    min.features = 200       # Cell must have at least 200 genes
)

porcn_ko <- CreateSeuratObject(
    counts      = porcn_ko, 
    project     = 'porcn_ko',
    min.cells   = 3,
    min.features = 200
)

# Add sample ID metadata
porcn_wt$ID <- "WT"
porcn_ko$ID <- "KO"

cat("✓ Seurat objects created\n")
cat("  WT:", ncol(porcn_wt), "cells\n")
cat("  KO:", ncol(porcn_ko), "cells\n")

################################################################################
# 3. Merge Objects
################################################################################

cat("\n=== Step 3: Merging Datasets ===\n")

porcn.combined <- merge(porcn_wt, y = porcn_ko)

# Verify merge
cat("✓ Merged object created:", ncol(porcn.combined), "total cells\n")
cat("  Sample distribution:\n")
print(table(porcn.combined$orig.ident))

# Set active identity
Idents(porcn.combined) <- 'ID'

################################################################################
# 4. Calculate QC Metrics
################################################################################

cat("\n=== Step 4: Calculating QC Metrics ===\n")

# Calculate mitochondrial percentage
# Note: Mouse mitochondrial genes start with "mt-" (lowercase)
porcn.combined[["percent.mt"]] <- PercentageFeatureSet(
    porcn.combined, 
    pattern = "^mt-"
)

cat("✓ Mitochondrial percentage calculated\n")

# Display summary statistics
cat("\nQC Metrics Summary:\n")
cat("  nFeature_RNA: ", summary(porcn.combined$nFeature_RNA), "\n")
cat("  nCount_RNA:   ", summary(porcn.combined$nCount_RNA), "\n")
cat("  percent.mt:   ", summary(porcn.combined$percent.mt), "\n")

################################################################################
# 5. Generate QC Visualizations
################################################################################

cat("\n=== Step 5: Generating QC Plots ===\n")

# Violin plots for QC metrics
vln_plot <- VlnPlot(
    porcn.combined,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0
) +
    plot_annotation(title = "QC Metrics Before Filtering")

ggsave(
    filename = "./results/01_QC/QC_violin_plots_before_filtering.png",
    plot     = vln_plot,
    width    = 12,
    height   = 5,
    dpi      = 300
)

# Scatter plots to examine relationships
plot1 <- FeatureScatter(
    porcn.combined,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
) + 
    ggtitle("UMI Count vs Mitochondrial %")

plot2 <- FeatureScatter(
    porcn.combined,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
) +
    ggtitle("UMI Count vs Gene Count")

scatter_plots <- plot1 + plot2

ggsave(
    filename = "./results/01_QC/QC_scatter_plots.png",
    plot     = scatter_plots,
    width    = 12,
    height   = 5,
    dpi      = 300
)

cat("✓ QC plots saved to ./results/01_QC/\n")

################################################################################
# 6. Apply Quality Filters
################################################################################

cat("\n=== Step 6: Applying Quality Filters ===\n")

# Count cells before filtering
n_before <- ncol(porcn.combined)

# Apply QC filters
porcn.combined <- subset(
    porcn.combined,
    subset = nFeature_RNA > QC_PARAMS$min_features & 
             nFeature_RNA < QC_PARAMS$max_features & 
             percent.mt < QC_PARAMS$max_mt_pct
)

# Count cells after filtering
n_after <- ncol(porcn.combined)
n_removed <- n_before - n_after
pct_removed <- round(100 * n_removed / n_before, 2)

cat("✓ Quality filtering complete\n")
cat("  Cells before filtering:", n_before, "\n")
cat("  Cells after filtering: ", n_after, "\n")
cat("  Cells removed:         ", n_removed, "(", pct_removed, "%)\n")

# Display filtered distribution
cat("\n  Filtered sample distribution:\n")
print(table(porcn.combined$ID))

################################################################################
# 7. Generate Post-Filtering QC Plots
################################################################################

cat("\n=== Step 7: Post-Filtering QC Plots ===\n")

vln_plot_filtered <- VlnPlot(
    porcn.combined,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0
) +
    plot_annotation(title = "QC Metrics After Filtering")

ggsave(
    filename = "./results/01_QC/QC_violin_plots_after_filtering.png",
    plot     = vln_plot_filtered,
    width    = 12,
    height   = 5,
    dpi      = 300
)

cat("✓ Post-filtering QC plots saved\n")

################################################################################
# 8. Save Processed Object
################################################################################

cat("\n=== Step 8: Saving Processed Object ===\n")

# Save the QC-filtered object
save(porcn.combined, file = "./data/porcn.combined_QC.RData")

cat("✓ QC-filtered object saved to ./data/porcn.combined_QC.RData\n")

################################################################################
# 9. Summary Report
################################################################################

cat("\n" %+% paste(rep("=", 80), collapse="") %+% "\n")
cat("QC SUMMARY REPORT\n")
cat(paste(rep("=", 80), collapse="") %+% "\n\n")

cat("Input Parameters:\n")
cat("  Min features per cell:      ", QC_PARAMS$min_features, "\n")
cat("  Max features per cell:      ", QC_PARAMS$max_features, "\n")
cat("  Max mitochondrial percent:  ", QC_PARAMS$max_mt_pct, "%\n\n")

cat("Filtering Results:\n")
cat("  Total cells (before):       ", n_before, "\n")
cat("  Total cells (after):        ", n_after, "\n")
cat("  Cells removed:              ", n_removed, "(", pct_removed, "%)\n\n")

cat("Final Object:\n")
cat("  Number of cells:            ", ncol(porcn.combined), "\n")
cat("  Number of genes:            ", nrow(porcn.combined), "\n")
cat("  WT cells:                   ", sum(porcn.combined$ID == "WT"), "\n")
cat("  KO cells:                   ", sum(porcn.combined$ID == "KO"), "\n\n")

cat("Output Files:\n")
cat("  - ./data/porcn.combined_QC.RData\n")
cat("  - ./results/01_QC/QC_violin_plots_before_filtering.png\n")
cat("  - ./results/01_QC/QC_violin_plots_after_filtering.png\n")
cat("  - ./results/01_QC/QC_scatter_plots.png\n\n")

cat(paste(rep("=", 80), collapse="") %+% "\n")
cat("✓ Step 01 Complete: Data Loading and Quality Control\n")
cat(paste(rep("=", 80), collapse="") %+% "\n\n")

# Clean up intermediate objects
rm(porcn_wt, porcn_ko)
gc()

cat("Next step: 02_integration/02_harmony_integration.R\n\n")
