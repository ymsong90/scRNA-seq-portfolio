################################################################################
# Single-Cell RNA-seq Analysis Pipeline
# Step 03: Clustering and Cell Type Annotation
# 
# Description:
#   - Find marker genes for each cluster
#   - Annotate cell types based on canonical markers
#   - Generate dotplots and feature plots for visualization
#
# Input:
#   - ./data/porcn.combined.harmony.RData (from Step 02)
#
# Output:
#   - Annotated Seurat object with cell type labels
#   - Marker gene tables
#   - Visualization plots
#
# Author: YMS
# Date: 2025
################################################################################

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)

# Create output directory
if (!dir.exists("./results/03_clustering")) {
    dir.create("./results/03_clustering", recursive = TRUE)
}

################################################################################
# 1. Load Integrated Data
################################################################################

load("./data/porcn.combined.harmony.RData")

################################################################################
# 2. Find Cluster Markers
################################################################################

# FindAllMarkers는 시간이 오래 걸릴 수 있음
# JoinLayers를 먼저 실행해야 Seurat v5에서 작동
porcn.combined.harmony <- JoinLayers(object = porcn.combined.harmony)

# NOTE: only.pos = TRUE는 upregulated markers만 찾음
# min.pct = 0.25는 최소 25% 세포에서 발현되어야 함
# logfc.threshold = 0.25는 최소 0.25 log fold-change
porcn.markers <- FindAllMarkers(
    porcn.combined.harmony,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
)

# Save marker genes
write.csv(
    porcn.markers,
    file = "./results/03_clustering/cluster_markers_all.csv",
    row.names = FALSE
)

# Display top markers per cluster
top_markers <- porcn.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

################################################################################
# 3. Visualize Canonical Cell Type Markers
################################################################################

# Define canonical markers for major cell types
canonical_markers <- list(
    # T cells & NK cells
    Tcell_NK = c('Ptprc', 'Cd3e', 'Cd3d', 'Cd4', 'Cd8a', 'Nkg7', 'Klrg1', 
                 'Sell', 'Prf1', 'Eomes', 'Trdc', 'Ccr7', 'Il7r', 'Cd28', 'Tbx21'),
    
    # Epithelial/Cancer cells
    Epithelial = c('Krt19', 'Cdh1', 'Epcam', 'Gli3', 'Tff1', 'Tff2', 'Tff3', 
                   'Mki67', 'Top2a', 'Ido1', 'Ido2', 'Msln'),
    
    # Myeloid cells
    Myeloid = c('Ptprc', 'Cd68', 'Adgre1', 'Itgam', 'Cd14', 'Mrc1', 'H2-Eb1', 
                'Batf3', 'S100a8', 'Ly6c2', 'Ly6g', 'Ccr7', 'Clec9a'),
    
    # Stromal cells
    Stromal = c('Col1a2', 'Pdpn', 'Pdgfa', 'Dcn', 'Cdh11', 'Apoe', 
                'Pecam1', 'Cdh5'),
    
    # B cells
    Bcell = c('Cd79a', 'Cd19', 'Ms4a1')
)

# Generate feature plots for each marker set
for (celltype in names(canonical_markers)) {
    features <- canonical_markers[[celltype]]
    
    p <- FeaturePlot(
        object   = porcn.combined.harmony,
        features = features,
        cols     = c("grey", "blue"),
        reduction = "umap",
        pt.size  = 0.1,
        ncol     = 4
    )
    
    ggsave(
        filename = paste0("./results/03_clustering/FeaturePlot_", celltype, ".png"),
        plot     = p,
        width    = 16,
        height   = ceiling(length(features) / 4) * 4,
        dpi      = 300
    )
}

################################################################################
# 4. Cell Type Annotation
################################################################################

# Store cluster assignments
porcn.combined.harmony[["UMAP_Clusters"]] <- Idents(object = porcn.combined.harmony)

# IMPORTANT: 실제 데이터의 cluster 번호에 맞춰 수정 필요
# Marker gene expression을 확인 후 적절한 cell type 할당
current.cluster.ids <- levels(Idents(porcn.combined.harmony))

# Cell type annotation mapping
new.cluster.ids <- c(
    "0"  = "Epi/Cancer cell",
    "1"  = "Epi/Cancer cell",
    "2"  = "Epi/Cancer cell",
    "3"  = "CTL",
    "4"  = "Epi/Cancer cell",
    "5"  = "Apoehi CAF",
    "6"  = "CD4 T cell",
    "7"  = "Mono/Mac",
    "8"  = "Dcnhi CAF",
    "9"  = "B cell",
    "10" = "Epi/Cancer cell",
    "11" = "NK cell",
    "12" = "Pdgfahi CAF",
    "13" = "Endothelial cell",
    "14" = "Epi/Cancer cell",
    "15" = "Mki67hiTop2ahi Cancer cell",
    "16" = "Batf3hi DC",
    "17" = "Neutrophil",
    "18" = "Epcamhi B cell",
    "19" = "Mki67hiTop2Ahi CTL",
    "20" = "Gli3hi Cancer cell",
    "21" = "Cdh11hi CAF",
    "22" = "Tffhi Cancer cell",
    "23" = "Mki67hiTop2Ahi B cell",
    "24" = "Mki67intTop2aint Cancer cell",
    "25" = "Mki67hiTop2AhiBatf3hi DC",
    "26" = "Idohi cancer cell"
)

# Apply annotations only for existing clusters
existing_clusters <- intersect(names(new.cluster.ids), current.cluster.ids)
annotations <- new.cluster.ids[existing_clusters]

porcn.combined.harmony <- RenameIdents(
    object = porcn.combined.harmony,
    !!!annotations
)

# Store annotations in metadata
porcn.combined.harmony[["NH_labels"]] <- Idents(porcn.combined.harmony)

# Display cell type distribution
table(porcn.combined.harmony$NH_labels)

################################################################################
# 5. Define Cell Type Order and Colors
################################################################################

# Define desired order for visualization
desired_order <- c(
    # Cancer cells (top)
    "Gli3hi Cancer cell", "Tffhi Cancer cell",
    "Mki67intTop2aint Cancer cell", "Mki67hiTop2ahi Cancer cell",
    "Idohi cancer cell", "Epi/Cancer cell",
    
    # T cells / NK cells
    "NK cell", "CD4 T cell", "Mki67hiTop2Ahi CTL", "CTL",
    
    # CAF
    "Pdgfahi CAF", "Dcnhi CAF", "Cdh11hi CAF", "Apoehi CAF",
    
    # Myeloid / Innate
    "Mki67hiTop2AhiBatf3hi DC", "Batf3hi DC", "Mono/Mac", "Neutrophil",
    
    # B cells
    "Mki67hiTop2Ahi B cell", "Epcamhi B cell", "B cell",
    
    # Endothelial (bottom)
    "Endothelial cell"
)

# Set factor levels
porcn.combined.harmony$NH_labels <- factor(
    porcn.combined.harmony$NH_labels,
    levels = desired_order
)

# Define color palette for cell types
celltype.cols <- c(
    # B-cell lineage
    "B cell"                    = "#808000",
    "Epcamhi B cell"            = "#B1D337",
    "Mki67hiTop2Ahi B cell"     = "#21A04A",
    
    # Neutrophil / DC / Macrophage
    "Neutrophil"                = "#778899",
    "Batf3hi DC"                = "#F17D97",
    "Mki67hiTop2AhiBatf3hi DC"  = "#A7213A",
    "Mono/Mac"                  = "#ED4169",
    
    # Cancer-cell lineage
    "Epi/Cancer cell"           = "#0B6A42",
    "Idohi cancer cell"         = "#CDE9C8",
    "Mki67intTop2aint Cancer cell" = "#0B6B6A",
    "Tffhi Cancer cell"         = "#22A1B1",
    "Mki67hiTop2ahi Cancer cell"= "#73CCD4",
    "Gli3hi Cancer cell"        = "#C2E5EB",
    
    # CAF lineage
    "Apoehi CAF"                = "#F65C4B",
    "Cdh11hi CAF"               = "#F68F80",
    "Dcnhi CAF"                 = "#EA6D22",
    "Pdgfahi CAF"               = "#F58927",
    
    # Endothelial
    "Endothelial cell"          = "#F0B817",
    
    # CTL / CD4 T
    "CTL"                       = "#C2E9F9",
    "Mki67hiTop2Ahi CTL"        = "#82D5F7",
    "CD4 T cell"                = "#0E7DC2",
    
    # NK
    "NK cell"                   = "#8F7BB7"
)

################################################################################
# 6. Generate Annotated UMAP
################################################################################

# Set identity to cell types
Idents(porcn.combined.harmony) <- "NH_labels"

# UMAP with cell type labels
umap_annotated <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    label     = TRUE,
    repel     = TRUE,
    label.size = 3,
    cols      = celltype.cols,
    pt.size   = 0.5
) +
    ggtitle("Annotated UMAP (All Cells)") +
    NoLegend()

ggsave(
    filename = "./results/03_clustering/UMAP_annotated.png",
    plot     = umap_annotated,
    width    = 10,
    height   = 8,
    dpi      = 300
)

# UMAP split by condition
umap_split <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    split.by  = "ID",
    label     = TRUE,
    repel     = TRUE,
    label.size = 2.5,
    cols      = celltype.cols,
    pt.size   = 0.3,
    ncol      = 2
) +
    ggtitle("Annotated UMAP (Split by Condition)")

ggsave(
    filename = "./results/03_clustering/UMAP_annotated_split.png",
    plot     = umap_split,
    width    = 16,
    height   = 7,
    dpi      = 300
)

# UMAP with legend
umap_with_legend <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    label     = FALSE,
    cols      = celltype.cols,
    pt.size   = 0.5
) +
    ggtitle("Annotated UMAP") +
    theme(legend.text = element_text(size = 10))

ggsave(
    filename = "./results/03_clustering/UMAP_annotated_with_legend.png",
    plot     = umap_with_legend,
    width    = 12,
    height   = 8,
    dpi      = 300
)

################################################################################
# 7. Generate DotPlot for Marker Genes
################################################################################

# Define marker genes for dotplot
total_features <- unique(c(
    # General markers
    'Ptprc', 'Krt19', 'Cdh1', 'Msln',
    # Cancer markers
    'Gli3', 'Tff1', 'Tff2', 'Tff3', 'Mki67', 'Top2a', 'Ido1', 'Ido2',
    # T cell markers
    'Nkg7', 'Klrg1', 'Cd3d', 'Cd3e', 'Cd4', 'Cd8a', 'Sell',
    'Ccr7', 'Il7r', 'Cd28', 'Il2ra', 'Il1rl1', 'Gata3', 'Trdc',
    # Stromal markers
    'Col1a2', 'Pdpn', 'Pdgfa', 'Dcn', 'Cdh11', 'Apoe',
    # Myeloid markers
    'Cd68', 'Adgre1', 'Itgam', 'Cd14', 'Mrc1', 'H2-Eb1', 'Batf3', 'S100a8',
    # B cell markers
    'Cd79a', 'Cd19', 'Ms4a1', 'Ly6c2', 'Ly6g',
    # Endothelial markers
    'Try4', 'Pecam1', 'Epcam', 'Cdh5'
))

dotplot <- DotPlot(
    object    = porcn.combined.harmony,
    features  = total_features,
    cols      = "RdYlBu",
    dot.scale = 10
) +
    RotatedAxis() +
    theme_bw(base_size = 12) +
    theme(
        axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                    colour = "black", size = 10),
        axis.text.y  = element_text(colour = "black", size = 10),
        axis.ticks   = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    ) +
    labs(x = NULL, y = NULL)

ggsave(
    filename = "./results/03_clustering/DotPlot_markers.png",
    plot     = dotplot,
    width    = 18,
    height   = 10,
    dpi      = 300
)

# High-resolution TIFF version
ggsave(
    filename = "./results/03_clustering/DotPlot_markers.tiff",
    plot     = dotplot,
    width    = 20,
    height   = 10,
    dpi      = 600,
    compression = "lzw"
)

################################################################################
# 8. Save Annotated Object
################################################################################

save(porcn.combined.harmony, file = "./data/porcn.combined.harmony_annotated.RData")

gc()

# NOTE: 다음 단계는 04_visualization/04_cell_proportion_analysis.R
