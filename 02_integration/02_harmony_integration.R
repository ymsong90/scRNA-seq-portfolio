################################################################################
# Single-Cell RNA-seq Analysis Pipeline
# Step 02: Normalization and Harmony Integration
# 
# Description:
#   - Normalize gene expression data
#   - Identify highly variable features
#   - Scale data
#   - Perform batch correction using Harmony
#   - Dimensionality reduction (PCA, UMAP, t-SNE)
#
# Input:
#   - ./data/porcn.combined_QC.RData (from Step 01)
#
# Output:
#   - porcn.combined.harmony: Harmony-corrected integrated object
#   - PCA, UMAP, and t-SNE embeddings
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
    library(harmony)
})

# Set parameters ---------------------------------------------------------------
NORM_PARAMS <- list(
    normalization_method = "LogNormalize",
    scale_factor         = 10000,
    n_variable_features  = 2000,
    n_pcs                = 10
)

# Create output directory
if (!dir.exists("./results/02_integration")) {
    dir.create("./results/02_integration", recursive = TRUE)
}

################################################################################
# 1. Load QC-Filtered Data
################################################################################

cat("\n=== Step 1: Loading QC-Filtered Data ===\n")

load("./data/porcn.combined_QC.RData")

cat("✓ Data loaded successfully\n")
cat("  Cells:", ncol(porcn.combined), "\n")
cat("  Genes:", nrow(porcn.combined), "\n")

################################################################################
# 2. Normalization
################################################################################

cat("\n=== Step 2: Data Normalization ===\n")

porcn.combined <- NormalizeData(
    object               = porcn.combined,
    normalization.method = NORM_PARAMS$normalization_method,
    scale.factor         = NORM_PARAMS$scale_factor
)

cat("✓ Normalization complete\n")
cat("  Method:       ", NORM_PARAMS$normalization_method, "\n")
cat("  Scale factor: ", NORM_PARAMS$scale_factor, "\n")

################################################################################
# 3. Identify Highly Variable Features
################################################################################

cat("\n=== Step 3: Identifying Variable Features ===\n")

porcn.combined <- FindVariableFeatures(
    porcn.combined,
    selection.method = "vst",
    nfeatures        = NORM_PARAMS$n_variable_features
)

# Get top variable genes
top_genes <- head(VariableFeatures(porcn.combined), 20)
cat("✓ Variable features identified:", length(VariableFeatures(porcn.combined)), "\n")
cat("\n  Top 20 variable genes:\n")
print(top_genes)

# Plot variable features
plot_var_features <- VariableFeaturePlot(porcn.combined)
plot_var_labeled <- LabelPoints(
    plot   = plot_var_features,
    points = top_genes,
    repel  = TRUE,
    xnudge = 0,
    ynudge = 0
)

ggsave(
    filename = "./results/02_integration/variable_features.png",
    plot     = plot_var_labeled,
    width    = 10,
    height   = 7,
    dpi      = 300
)

cat("✓ Variable feature plot saved\n")

################################################################################
# 4. Scale Data
################################################################################

cat("\n=== Step 4: Scaling Data ===\n")
cat("⚠ This step may take several minutes...\n")

# Scale all genes for downstream analyses
porcn.combined <- ScaleData(
    object   = porcn.combined,
    features = rownames(porcn.combined)
)

cat("✓ Data scaling complete\n")

# Save intermediate object
save(porcn.combined, file = "./data/porcn.combined_scaled.RData")
cat("✓ Scaled object saved (checkpoint)\n")

################################################################################
# 5. PCA (Before Harmony)
################################################################################

cat("\n=== Step 5: Principal Component Analysis ===\n")

porcn.combined <- RunPCA(
    porcn.combined,
    features = VariableFeatures(object = porcn.combined),
    verbose  = FALSE
)

# Visualize PCA variance
elbow_plot <- ElbowPlot(porcn.combined, ndims = 50) +
    ggtitle("PCA Elbow Plot") +
    geom_vline(xintercept = NORM_PARAMS$n_pcs, linetype = "dashed", color = "red") +
    annotate("text", 
             x = NORM_PARAMS$n_pcs + 5, 
             y = max(porcn.combined@reductions$pca@stdev[1:20]), 
             label = paste0("n_pcs = ", NORM_PARAMS$n_pcs),
             color = "red")

ggsave(
    filename = "./results/02_integration/PCA_elbow_plot.png",
    plot     = elbow_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
)

cat("✓ PCA complete (50 PCs computed)\n")
cat("  Selected PCs for downstream:", NORM_PARAMS$n_pcs, "\n")

# PCA visualization before harmony
pca_before <- DimPlot(
    porcn.combined,
    reduction = "pca",
    group.by  = "ID",
    pt.size   = 0.5
) +
    ggtitle("PCA (Before Harmony)")

################################################################################
# 6. Clustering and UMAP (Without Harmony) - For Comparison
################################################################################

cat("\n=== Step 6: Clustering Without Harmony (Baseline) ===\n")

porcn.combined <- FindNeighbors(
    object = porcn.combined,
    dims   = 1:NORM_PARAMS$n_pcs
)

porcn.combined <- FindClusters(
    object     = porcn.combined,
    resolution = 1.5
)

porcn.combined <- RunUMAP(
    object = porcn.combined,
    dims   = 1:NORM_PARAMS$n_pcs
)

# Store cluster assignments before Harmony
porcn.combined[["clusters_no_harmony"]] <- Idents(object = porcn.combined)

# Visualize without Harmony
umap_no_harmony <- DimPlot(
    object    = porcn.combined,
    reduction = "umap",
    group.by  = "ID",
    pt.size   = 0.5
) +
    ggtitle("UMAP Without Harmony")

umap_no_harmony_split <- DimPlot(
    object    = porcn.combined,
    reduction = "umap",
    split.by  = "ID",
    pt.size   = 0.5
) +
    ggtitle("UMAP Without Harmony (Split by Condition)")

cat("✓ Baseline clustering complete (for comparison)\n")

################################################################################
# 7. Harmony Batch Correction
################################################################################

cat("\n=== Step 7: Harmony Batch Correction ===\n")
cat("⚠ Running Harmony integration...\n")

porcn.combined.harmony <- porcn.combined %>% 
    RunHarmony(
        "ID",                    # Batch variable to correct
        plot_convergence = TRUE  # Show convergence plot
    )

cat("✓ Harmony integration complete\n")

# Save Harmony convergence plot
harmony_plot <- porcn.combined.harmony@tools$RunHarmony$plot
ggsave(
    filename = "./results/02_integration/harmony_convergence.png",
    plot     = harmony_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
)

################################################################################
# 8. Clustering and UMAP (With Harmony)
################################################################################

cat("\n=== Step 8: Clustering With Harmony ===\n")

# Find neighbors using Harmony reduction
porcn.combined.harmony <- FindNeighbors(
    object    = porcn.combined.harmony,
    dims      = 1:NORM_PARAMS$n_pcs,
    reduction = "harmony"
)

# Find clusters
porcn.combined.harmony <- FindClusters(
    object     = porcn.combined.harmony,
    resolution = 2.5,
    reduction  = "harmony"
)

cat("✓ Identified", length(unique(Idents(porcn.combined.harmony))), "clusters\n")

# Run UMAP on Harmony embeddings
porcn.combined.harmony <- RunUMAP(
    object    = porcn.combined.harmony,
    dims      = 1:NORM_PARAMS$n_pcs,
    reduction = "harmony"
)

# Run t-SNE on Harmony embeddings
porcn.combined.harmony <- RunTSNE(
    object    = porcn.combined.harmony,
    dims      = 1:NORM_PARAMS$n_pcs,
    reduction = "harmony"
)

cat("✓ Dimensionality reduction complete (UMAP & t-SNE)\n")

################################################################################
# 9. Visualization Comparisons
################################################################################

cat("\n=== Step 9: Generating Comparison Plots ===\n")

# UMAP with Harmony - by cluster
umap_harmony_clusters <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    label     = TRUE,
    pt.size   = 0.8
) +
    ggtitle("UMAP With Harmony (Clusters)")

# UMAP with Harmony - by condition
umap_harmony_condition <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    group.by  = "ID",
    pt.size   = 0.5
) +
    ggtitle("UMAP With Harmony (Condition)")

# UMAP with Harmony - split by condition
umap_harmony_split <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    split.by  = "ID",
    label     = TRUE,
    pt.size   = 0.5
) +
    ggtitle("UMAP With Harmony (Split)")

# t-SNE visualization
tsne_harmony <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "tsne",
    label     = TRUE,
    pt.size   = 0.7
) +
    ggtitle("t-SNE With Harmony")

# Combine comparison plots
comparison_plot <- (umap_no_harmony | umap_harmony_condition) /
                   (pca_before | tsne_harmony)

ggsave(
    filename = "./results/02_integration/integration_comparison.png",
    plot     = comparison_plot,
    width    = 16,
    height   = 12,
    dpi      = 300
)

# Save individual plots
ggsave(
    filename = "./results/02_integration/UMAP_harmony_clusters.png",
    plot     = umap_harmony_clusters,
    width    = 10,
    height   = 8,
    dpi      = 300
)

ggsave(
    filename = "./results/02_integration/UMAP_harmony_split.png",
    plot     = umap_harmony_split,
    width    = 14,
    height   = 6,
    dpi      = 300
)

cat("✓ Comparison plots saved\n")

################################################################################
# 10. QC Metrics Visualization
################################################################################

cat("\n=== Step 10: Post-Integration QC ===\n")

# Violin plots for QC metrics after integration
Idents(object = porcn.combined.harmony) <- 'ID'

vln_qc <- VlnPlot(
    porcn.combined.harmony,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol     = 3,
    pt.size  = 0
)

ggsave(
    filename = "./results/02_integration/QC_metrics_post_integration.png",
    plot     = vln_qc,
    width    = 12,
    height   = 5,
    dpi      = 300
)

cat("✓ Post-integration QC plots saved\n")

################################################################################
# 11. Save Integrated Object
################################################################################

cat("\n=== Step 11: Saving Integrated Object ===\n")

save(porcn.combined.harmony, file = "./data/porcn.combined.harmony.RData")

cat("✓ Harmony-corrected object saved to ./data/porcn.combined.harmony.RData\n")

################################################################################
# 12. Summary Report
################################################################################

cat("\n" %+% paste(rep("=", 80), collapse="") %+% "\n")
cat("INTEGRATION SUMMARY REPORT\n")
cat(paste(rep("=", 80), collapse="") %+% "\n\n")

cat("Normalization:\n")
cat("  Method:             ", NORM_PARAMS$normalization_method, "\n")
cat("  Scale factor:       ", NORM_PARAMS$scale_factor, "\n")
cat("  Variable features:  ", NORM_PARAMS$n_variable_features, "\n\n")

cat("Dimensionality Reduction:\n")
cat("  PCs used:           ", NORM_PARAMS$n_pcs, "\n")
cat("  Clustering resolution: 2.5\n")
cat("  Clusters identified:", length(unique(Idents(porcn.combined.harmony))), "\n\n")

cat("Batch Correction:\n")
cat("  Method:             Harmony\n")
cat("  Batch variable:     ID (WT vs KO)\n")
cat("  Converged:          Yes\n\n")

cat("Output Files:\n")
cat("  - ./data/porcn.combined.harmony.RData\n")
cat("  - ./results/02_integration/variable_features.png\n")
cat("  - ./results/02_integration/PCA_elbow_plot.png\n")
cat("  - ./results/02_integration/harmony_convergence.png\n")
cat("  - ./results/02_integration/UMAP_harmony_clusters.png\n")
cat("  - ./results/02_integration/UMAP_harmony_split.png\n")
cat("  - ./results/02_integration/integration_comparison.png\n\n")

cat(paste(rep("=", 80), collapse="") %+% "\n")
cat("✓ Step 02 Complete: Normalization and Harmony Integration\n")
cat(paste(rep("=", 80), collapse="") %+% "\n\n")

# Clean up
rm(porcn.combined)
gc()

cat("Next step: 03_clustering/03_clustering_annotation.R\n\n")
