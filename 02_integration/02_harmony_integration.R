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

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

# Set parameters
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

load("./data/porcn.combined_QC.RData")

################################################################################
# 2. Normalization
################################################################################

# NOTE: LogNormalize는 각 cell의 총 카운트를 10,000으로 정규화 후 log transformation
porcn.combined <- NormalizeData(
    object               = porcn.combined,
    normalization.method = NORM_PARAMS$normalization_method,
    scale.factor         = NORM_PARAMS$scale_factor
)

################################################################################
# 3. Identify Highly Variable Features
################################################################################

# NOTE: vst method는 variance-to-mean ratio를 이용한 방법
# 생물학적 변이가 큰 유전자를 선별
porcn.combined <- FindVariableFeatures(
    porcn.combined,
    selection.method = "vst",
    nfeatures        = NORM_PARAMS$n_variable_features
)

# Get top variable genes
top_genes <- head(VariableFeatures(porcn.combined), 20)

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

################################################################################
# 4. Scale Data
################################################################################

# IMPORTANT: 모든 유전자를 scaling (평균 0, 분산 1로 표준화)
# Downstream analysis (PCA, clustering 등)를 위해 필수
porcn.combined <- ScaleData(
    object   = porcn.combined,
    features = rownames(porcn.combined)
)

# Save intermediate object (checkpoint)
save(porcn.combined, file = "./data/porcn.combined_scaled.RData")

################################################################################
# 5. PCA (Before Harmony)
################################################################################

# NOTE: PCA는 dimensionality reduction의 첫 단계
# Variable features만 사용하여 계산
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

# NOTE: Harmony 적용 전 baseline 결과 확인용
# Batch effect가 있는지 확인하기 위한 비교군
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

################################################################################
# 7. Harmony Batch Correction
################################################################################

# IMPORTANT: Harmony는 iterative clustering 기반 batch correction
# "ID" 변수 (WT vs KO)에 대한 batch effect를 보정
# plot_convergence = TRUE로 수렴 과정 확인 가능
porcn.combined.harmony <- porcn.combined %>% 
    RunHarmony(
        "ID",
        plot_convergence = TRUE
    )

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

# NOTE: Harmony-corrected embedding을 사용
# reduction = "harmony"로 지정
porcn.combined.harmony <- FindNeighbors(
    object    = porcn.combined.harmony,
    dims      = 1:NORM_PARAMS$n_pcs,
    reduction = "harmony"
)

# 2.5는 상대적으로 높은 값 (fine-grained clustering)
porcn.combined.harmony <- FindClusters(
    object     = porcn.combined.harmony,
    resolution = 2.5,
    reduction  = "harmony"
)

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

################################################################################
# 9. Visualization Comparisons
################################################################################

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

################################################################################
# 10. QC Metrics Visualization
################################################################################

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

################################################################################
# 11. Save Integrated Object
################################################################################

save(porcn.combined.harmony, file = "./data/porcn.combined.harmony.RData")

# Clean up
rm(porcn.combined)
gc()

# NOTE: 다음 단계는 03_clustering/03_clustering_annotation.R
