################################################################################
# IMC Analysis Pipeline
# Step 02: Batch Correction and Clustering
# 
# Description:
#   - Run PCA on marker expression
#   - Perform Harmony batch correction by sample_id
#   - SOM clustering (100 codes)
#   - ConsensusClusterPlus metaclustering (K=35)
#   - Generate cluster-level heatmaps
#   - UMAP on Harmony-corrected embeddings
#
# Input:
#   - ./data/spe_preprocessed.rds (from Step 01)
#
# Output:
#   - ./data/spe_batch_corrected.rds
#   - Cluster heatmaps and UMAP plots
#
# Author: YMS
# Date: 2024-11-17
################################################################################

# Load required packages
library(harmony)
library(BiocSingular)
library(scater)
library(dittoSeq)
library(viridis)
library(cowplot)
library(SpatialExperiment)
library(bluster)
library(ConsensusClusterPlus)
library(scran)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(dplyr)

# Set parameters
BATCH_CORRECTION_PARAMS <- list(
    n_pcs              = 30,           # Number of PCs for Harmony
    som_codes          = 100,          # Number of SOM codes
    n_clusters         = 35,           # Final number of clusters
    ccp_reps           = 100,          # ConsensusClusterPlus repetitions
    ccp_max_k          = 60,           # Maximum K for ConsensusClusterPlus
    seed               = 900223
)

# Create output directory
if (!dir.exists("./figure/02_batch_correction")) {
    dir.create("./figure/02_batch_correction", recursive = TRUE)
}

################################################################################
# 1. Load Preprocessed Data
################################################################################

spe <- readRDS("./data/spe_preprocessed.rds")

# Ensure use_channel is defined
if (!"use_channel" %in% colnames(rowData(spe))) {
    exclude_markers <- c("DNA1", "DNA2")
    rowData(spe)$use_channel <- !grepl(
        paste(exclude_markers, collapse = "|"),
        rownames(spe),
        ignore.case = TRUE
    )
}

################################################################################
# 2. PCA on Marker Expression
################################################################################

# NOTE: PCA is performed on transformed expression values
# Only use markers defined in use_channel
set.seed(BATCH_CORRECTION_PARAMS$seed)

spe <- runPCA(
    spe,
    subset_row = rowData(spe)$use_channel,
    exprs_values = "exprs",
    ncomponents = BATCH_CORRECTION_PARAMS$n_pcs,
    BSPARAM = ExactParam()
)

################################################################################
# 3. Harmony Batch Correction
################################################################################

# IMPORTANT: Harmony corrects batch effects while preserving biological variation
# Correction is performed on PCA space
# group.by.vars should be the batch variable (typically sample_id)
set.seed(BATCH_CORRECTION_PARAMS$seed)

out <- RunHarmony(
    spe, 
    group.by.vars = "sample_id"
)

# Extract Harmony embeddings
# NOTE: Harmony package may use "HARMONY" or "harmony" as reduction name
harm_name <- if ("HARMONY" %in% reducedDimNames(out)) "HARMONY" else "harmony"
reducedDim(spe, "harmony") <- reducedDim(out, harm_name)

################################################################################
# 4. UMAP on Harmony Embeddings
################################################################################

set.seed(BATCH_CORRECTION_PARAMS$seed)

spe <- runUMAP(
    spe, 
    dimred = "harmony", 
    name = "UMAP_harmonyCorrected"
)

# Optional: t-SNE on Harmony embeddings
# spe <- runTSNE(
#     spe,
#     dimred = "harmony",
#     name = "TSNE_harmonyCorrected"
# )

################################################################################
# 5. Visualization: Before vs After Correction
################################################################################

p_before <- dittoDimPlot(
    spe, 
    var = "group", 
    reduction.use = "UMAP", 
    size = 0.2
) +
    ggtitle("UMAP Before Harmony") +
    scale_color_manual(values = metadata(spe)$color_vectors$group)

p_after <- dittoDimPlot(
    spe, 
    var = "group", 
    reduction.use = "UMAP_harmonyCorrected", 
    size = 0.2
) +
    ggtitle("UMAP After Harmony") +
    scale_color_manual(values = metadata(spe)$color_vectors$group)

ggsave(
    filename = "./figure/02_batch_correction/01_harmony_comparison.tiff",
    plot = plot_grid(p_before, p_after, ncol = 2),
    width = 14,
    height = 6,
    dpi = 300,
    compression = "lzw"
)

################################################################################
# 6. SOM Clustering
################################################################################

# Extract Harmony embedding matrix
mat <- reducedDim(spe, "harmony")

# NOTE: SOM creates a grid of prototype vectors (codes)
# This provides an initial grouping of cells
som.out <- clusterRows(
    mat, 
    SomParam(BATCH_CORRECTION_PARAMS$som_codes), 
    full = TRUE
)

################################################################################
# 7. ConsensusClusterPlus Metaclustering
################################################################################

# IMPORTANT: ConsensusClusterPlus performs hierarchical clustering
# on SOM codes with consensus approach for stability
set.seed(BATCH_CORRECTION_PARAMS$seed)

ccp <- ConsensusClusterPlus(
    t(som.out$objects$som$codes[[1]]),  # Transpose: codes as rows
    maxK = BATCH_CORRECTION_PARAMS$ccp_max_k,
    reps = BATCH_CORRECTION_PARAMS$ccp_reps,
    distance = "euclidean",
    seed = BATCH_CORRECTION_PARAMS$seed,
    plot = NULL
)

# Select desired K
K <- BATCH_CORRECTION_PARAMS$n_clusters
som.cluster <- ccp[[K]][["consensusClass"]][som.out$clusters]
spe$som_clusters_harmonyCorrected <- factor(som.cluster)

################################################################################
# 8. Visualization: Clusters on UMAP
################################################################################

p_clusters <- dittoDimPlot(
    spe, 
    var = "som_clusters_harmonyCorrected",
    reduction.use = "UMAP_harmonyCorrected",
    size = 0.2, 
    do.label = TRUE
) +
    ggtitle(sprintf("SOM + ConsensusClusterPlus (K=%d) on Harmony UMAP", K))

ggsave(
    filename = "./figure/02_batch_correction/02_clusters_UMAP.tiff",
    plot = p_clusters,
    width = 10,
    height = 8,
    dpi = 300,
    compression = "lzw"
)

################################################################################
# 9. Define Cluster Color Palette
################################################################################

# Create color palette for all clusters
make_cluster_colors <- function(n) {
    pals <- c("Set3", "Paired", "Set1", "Set2", "Accent", "Dark2", "Pastel1", "Pastel2")
    cols <- unlist(lapply(pals, function(p) {
        k <- brewer.pal.info[p, "maxcolors"]
        brewer.pal(as.integer(k), p)
    }))
    cols <- unique(cols)
    if (length(cols) < n) {
        cols <- c(cols, hue_pal()(n - length(cols)))
    }
    cols[seq_len(n)]
}

# Create cluster celltype column
spe$cluster_celltype <- factor(as.integer(spe$som_clusters_harmonyCorrected))
spe$outline_celltype <- spe$cluster_celltype

levs <- levels(spe$cluster_celltype)
n_lev <- length(levs)

paletteN <- make_cluster_colors(max(BATCH_CORRECTION_PARAMS$n_clusters, n_lev))

cluster_cols <- setNames(paletteN[seq_len(n_lev)], levs)
outline_cols <- setNames(rep("#000000", n_lev), levs)

# Store in metadata
metadata(spe)$color_vectors$cluster_celltype <- cluster_cols
metadata(spe)$color_vectors$outline_celltype <- outline_cols

################################################################################
# 10. Define Marker Classes
################################################################################

# Type markers: cell identity markers
type_markers <- c(
    "SMA", "CD45", "Vimentin", "CD31", "MHC-II", "CD11c", "CD68", "CD11b",
    "PanCK", "CD4", "B220", "CD8", "Trem2", "EpCAM", "Col1A1",
    "CD3", "CCR2", "Arg1", "Hexb", "CX3CR1", "CD206", "PDGFRa"
)

# State markers: functional markers
state_markers <- c(
    "b-catenin", "non-p-b-catenin", "p-c-Jun", "Wnt5a", "p-JNK",
    "GranzymeB", "MPO", "IL-1b", "iNOS",
    "VEGFA", "Ki67", "CleavCaspase3"
)

# Assign marker classes
mc <- rep("other", nrow(spe))
names(mc) <- rownames(spe)
mc[rownames(spe) %in% type_markers] <- "type"
mc[rownames(spe) %in% state_markers] <- "state"
rowData(spe)$marker_class <- mc

################################################################################
# 11. Cluster-Level Mean Heatmap
################################################################################

# Calculate mean expression per cluster
keep_rows <- rownames(spe)[rowData(spe)$marker_class %in% c("type", "state")]

celltype_mean <- aggregateAcrossCells(
    as(spe, "SingleCellExperiment"),
    ids = spe$cluster_celltype,
    statistics = "mean",
    use.assay.type = "exprs",
    subset.row = keep_rows
)

# Heatmap: No scaling
tiff(
    filename = "./figure/02_batch_correction/03_cluster_heatmap_no_scale.tiff",
    width = 13,
    height = 9,
    units = "in",
    res = 600,
    compression = "lzw"
)
dittoHeatmap(
    celltype_mean,
    assay = "exprs",
    cluster_cols = TRUE,
    scale = "none",
    heatmap.colors = viridis(100),
    annot.by = c("cluster_celltype", "ncells"),
    annotation_colors = list(
        cluster_celltype = metadata(spe)$color_vectors$cluster_celltype,
        ncells = plasma(100)
    )
)
dev.off()

# Heatmap: Scaled to max
tiff(
    filename = "./figure/02_batch_correction/04_cluster_heatmap_scaled_to_max.tiff",
    width = 13,
    height = 9,
    units = "in",
    res = 600,
    compression = "lzw"
)
dittoHeatmap(
    celltype_mean,
    assay = "exprs",
    cluster_cols = TRUE,
    scaled.to.max = TRUE,
    heatmap.colors.max.scaled = viridis(100, option = "C"),
    annot.by = c("cluster_celltype", "ncells"),
    annotation_colors = list(
        cluster_celltype = metadata(spe)$color_vectors$cluster_celltype,
        ncells = plasma(100)
    ),
    border_color = "grey80"
)
dev.off()

# Heatmap: Z-score scaling
tiff(
    filename = "./figure/02_batch_correction/05_cluster_heatmap_zscore.tiff",
    width = 10,
    height = 13,
    units = "in",
    res = 300,
    compression = "lzw"
)
dittoHeatmap(
    celltype_mean,
    assay = "exprs",
    cluster_cols = TRUE,
    scale = "row",
    breaks = seq(-4, 4, length.out = 101),
    legend_breaks = c(-4, -2, 0, 2, 4),
    annot.by = c("cluster_celltype", "ncells"),
    annotation_colors = list(
        cluster_celltype = metadata(spe)$color_vectors$cluster_celltype,
        ncells = plasma(100)
    )
)
dev.off()

################################################################################
# 12. Final UMAP with Cluster Colors
################################################################################

p_final <- dittoDimPlot(
    spe, 
    var = "cluster_celltype",
    reduction.use = "UMAP_harmonyCorrected",
    size = 0.5,
    do.label = FALSE
) +
    scale_color_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank()
    )

ggsave(
    filename = "./figure/02_batch_correction/06_UMAP_final_clusters.tiff",
    plot = p_final,
    width = 9,
    height = 7,
    dpi = 400,
    compression = "lzw"
)

################################################################################
# 13. Save Object
################################################################################

saveRDS(spe, "./data/spe_batch_corrected.rds")

# Clean up
rm(mat, som.out, ccp, celltype_mean, out)
gc()

# NOTE: Next step is 03_clustering_and_annotation.R
