################################################################################
# CellChat Analysis: Cell-Cell Communication Network
#
# Dataset: Human Pancreatic Cancer scRNA-seq
# Purpose: Analyze WNT signaling pathway-mediated cell-cell interactions
#          Focus on Myeloid → Epithelial cell communication
#
# Author: YMS
# Date: 2025-07
################################################################################

library(Seurat)
library(CellChat)
library(ggplot2)
library(stringr)
library(dplyr)

# Output directory
output_dir <- "./results/CellChat"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

################################################################################
# 1. Prepare Input Data
################################################################################

# Extract normalized expression data
data.input <- GetAssayData(seu[["RNA"]], layer = "data")

# NOTE: Human gene symbols need formatting
# - Remove Ensembl version numbers (e.g., ".1", ".2")
# - Convert to uppercase for CellChatDB compatibility
rownames(data.input) <- toupper(str_remove(rownames(data.input), "\\.[0-9]+$"))

# Prepare metadata
meta <- seu@meta.data
meta$samples <- factor("sample1")  # Required by CellChat

################################################################################
# 2. Create CellChat Object
################################################################################

# Create CellChat object with cell type annotation
cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "Celltype"
)

# Use full human ligand-receptor database
cellchat@DB <- CellChatDB.human

################################################################################
# 3. Standard CellChat Pipeline
################################################################################

# Subset data for analysis
cellchat <- subsetData(cellchat)

# Identify overexpressed ligands/receptors
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability
cellchat <- computeCommunProb(cellchat)

# Filter low-confidence interactions
# NOTE: min.cells threshold affects sensitivity
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Aggregate pathway-level communication
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Save processed CellChat object
saveRDS(cellchat, file.path(output_dir, "cellchat_processed.rds"))

################################################################################
# 4. WNT Signaling Network - All Cell Types
################################################################################

# Circle plot showing all WNT-mediated interactions
p_wnt_all <- netVisual_aggregate(
    cellchat,
    signaling = "WNT",
    layout = "circle",
    vertex.weight = "auto",    # Node size by cell abundance
    weight.scale = TRUE,       # Edge width by interaction strength
    label.edge = FALSE
)

ggsave(
    file.path(output_dir, "WNT_All_Celltypes_circle.tiff"),
    plot = p_wnt_all,
    width = 7, height = 7, dpi = 300
)

################################################################################
# 5. WNT Signaling Network - Myeloid to Epithelial
################################################################################

# Focus on Myeloid → Epithelial cell WNT signaling
# RATIONALE: Myeloid cells may regulate cancer cell behavior via WNT
p_wnt_me <- netVisual_aggregate(
    cellchat,
    signaling = "WNT",
    sources.use = "Myeloid",
    targets.use = "Epithelial cell",
    layout = "chord"
)

ggsave(
    file.path(output_dir, "WNT_Myeloid_to_Epithelial_chord.pdf"),
    plot = p_wnt_me,
    width = 8, height = 6
)

################################################################################
# 6. Export Communication Results
################################################################################

# Extract all WNT interactions
wnt_all <- subsetCommunication(
    cellchat,
    signaling = "WNT",
    slot.name = "net"
)

write.csv(
    wnt_all,
    file.path(output_dir, "WNT_All_Celltypes_interactions.csv"),
    row.names = FALSE
)

# Extract Myeloid → Epithelial WNT interactions
wnt_myeloid_epithelial <- wnt_all %>%
    filter(source == "Myeloid" & target == "Epithelial cell")

write.csv(
    wnt_myeloid_epithelial,
    file.path(output_dir, "WNT_Myeloid_to_Epithelial_interactions.csv"),
    row.names = FALSE
)

################################################################################
# End of Analysis
################################################################################
