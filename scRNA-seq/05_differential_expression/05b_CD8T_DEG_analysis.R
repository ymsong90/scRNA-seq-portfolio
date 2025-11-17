################################################################################
# Step 05b: CD8 T Cell Differential Expression Analysis
#
# Purpose: Analyze CD8 T cell transcriptional changes (KO vs WT)
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./data/porcn.combined.harmony_annotated.RData
#
# Output:
#   - ./results/05_DEG/cd8/*.csv (DEG tables)
#   - ./results/05_DEG/cd8/volcano/*.png
#
# Author: YMS
# Date: 2025
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Parameters
DEG_PARAMS <- list(
    logfc_threshold = 0.25,
    min_pct = 0.10,
    padj_cutoff = 0.05
)

# Create output directories
dir.create("./results/05_DEG/cd8", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/05_DEG/cd8/volcano", recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. Load Annotated Data
################################################################################

load("./data/porcn.combined.harmony_annotated.RData")

################################################################################
# 2. Extract CD8 T Cells
################################################################################

# NOTE: CD8 T cell includes CTL and Mki67hi CTL subtypes
cd8_cells <- subset(
    porcn.combined.harmony,
    subset = Complete_Labels == "CD8 T cell"
)

################################################################################
# 3. DEG Analysis: CD8 KO vs WT
################################################################################

# Prepare identity for comparison
cd8_cells$celltype_condition <- paste(
    cd8_cells$Complete_Labels,
    cd8_cells$ID,
    sep = "_"
)

Idents(cd8_cells) <- "celltype_condition"

# Find markers: KO vs WT within CD8 T cells
deg_cd8 <- FindMarkers(
    cd8_cells,
    ident.1 = "CD8 T cell_KO",
    ident.2 = "CD8 T cell_WT",
    min.pct = DEG_PARAMS$min_pct,
    logfc.threshold = DEG_PARAMS$logfc_threshold,
    test.use = "wilcox"
)

deg_cd8 <- deg_cd8 %>%
    rownames_to_column("gene")

write.csv(
    deg_cd8,
    file = "./results/05_DEG/cd8/DEG_CD8_KO_vs_WT.csv",
    row.names = FALSE
)

################################################################################
# 4. DEG Analysis: CD8 vs Other Cell Types
################################################################################

# NOTE: This identifies CD8 T cell-specific markers
Idents(porcn.combined.harmony) <- "Complete_Labels"

deg_cd8_vs_rest <- FindMarkers(
    porcn.combined.harmony,
    ident.1 = "CD8 T cell",
    ident.2 = NULL,  # Compare against all other cells
    min.pct = DEG_PARAMS$min_pct,
    logfc.threshold = DEG_PARAMS$logfc_threshold,
    test.use = "wilcox"
)

deg_cd8_vs_rest <- deg_cd8_vs_rest %>%
    rownames_to_column("gene")

write.csv(
    deg_cd8_vs_rest,
    file = "./results/05_DEG/cd8/DEG_CD8_vs_Rest.csv",
    row.names = FALSE
)

################################################################################
# 5. Generate Volcano Plots
################################################################################

# Volcano plot function
plot_volcano <- function(deg_df, title,
                         lfc_cut = DEG_PARAMS$logfc_threshold,
                         padj_cut = DEG_PARAMS$padj_cutoff,
                         label_n = 20) {
    
    deg_df <- deg_df %>%
        mutate(
            log10P = -log10(p_val_adj + 1e-300),
            sig_flag = case_when(
                p_val_adj < padj_cut & avg_log2FC > lfc_cut ~ "Up",
                p_val_adj < padj_cut & avg_log2FC < -lfc_cut ~ "Down",
                TRUE ~ "NS"
            )
        )
    
    top_genes <- deg_df %>%
        filter(sig_flag != "NS") %>%
        arrange(p_val_adj) %>%
        group_by(sig_flag) %>%
        slice_head(n = label_n/2)
    
    ggplot(deg_df, aes(x = avg_log2FC, y = log10P, color = sig_flag)) +
        geom_point(size = 1.5, alpha = 0.6) +
        scale_color_manual(
            values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "grey50")
        ) +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), 
                   linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(padj_cut), 
                   linetype = "dashed", alpha = 0.5) +
        geom_text_repel(
            data = top_genes,
            aes(label = gene),
            size = 3,
            max.overlaps = 20,
            show.legend = FALSE
        ) +
        labs(
            title = title,
            x = "log2 Fold Change",
            y = "-log10(Adjusted P-value)"
        ) +
        theme_bw(base_size = 12)
}

# Volcano: CD8 KO vs WT
p1 <- plot_volcano(deg_cd8, "CD8 T cell: KO vs WT")
ggsave(
    "./results/05_DEG/cd8/volcano/Volcano_CD8_KO_vs_WT.png",
    plot = p1, width = 10, height = 8, dpi = 300
)

# Volcano: CD8 vs Rest
p2 <- plot_volcano(deg_cd8_vs_rest, "CD8 T cell vs Other Cell Types")
ggsave(
    "./results/05_DEG/cd8/volcano/Volcano_CD8_vs_Rest.png",
    plot = p2, width = 10, height = 8, dpi = 300
)

################################################################################
# 6. Filter Significant DEGs
################################################################################

# CD8 KO vs WT - upregulated
deg_cd8_up <- deg_cd8 %>%
    filter(p_val_adj < DEG_PARAMS$padj_cutoff & 
           avg_log2FC > DEG_PARAMS$logfc_threshold)

write.csv(deg_cd8_up, 
          "./results/05_DEG/cd8/DEG_CD8_upregulated_KO.csv",
          row.names = FALSE)

# CD8 KO vs WT - downregulated
deg_cd8_down <- deg_cd8 %>%
    filter(p_val_adj < DEG_PARAMS$padj_cutoff & 
           avg_log2FC < -DEG_PARAMS$logfc_threshold)

write.csv(deg_cd8_down,
          "./results/05_DEG/cd8/DEG_CD8_downregulated_KO.csv",
          row.names = FALSE)

# Clean up
rm(cd8_cells, deg_cd8, deg_cd8_vs_rest, deg_cd8_up, deg_cd8_down, p1, p2)
gc()

################################################################################
# End of Step 05b
################################################################################
