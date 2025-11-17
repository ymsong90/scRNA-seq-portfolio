################################################################################
# Step 06a: Myeloid Cell GO Enrichment Analysis
#
# Purpose: Functional enrichment analysis of myeloid DEGs
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./results/05_DEG/myeloid/DEG_upregulated.csv
#   - ./results/05_DEG/myeloid/DEG_downregulated.csv
#
# Output:
#   - ./results/06_GO/myeloid/*.csv (GO results)
#   - ./results/06_GO/myeloid/*.png (GO plots)
#
# Author: YMS
# Date: 2025
################################################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Create output directory
dir.create("./results/06_GO/myeloid", recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. Load DEG Results
################################################################################

deg_up <- read.csv("./results/05_DEG/myeloid/DEG_upregulated.csv")
deg_down <- read.csv("./results/05_DEG/myeloid/DEG_downregulated.csv")

################################################################################
# 2. Convert Gene Symbols to Entrez IDs
################################################################################

# NOTE: GO enrichment requires Entrez gene IDs
convert_to_entrez <- function(gene_symbols) {
    entrez <- bitr(
        gene_symbols,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Mm.eg.db
    )
    return(entrez$ENTREZID)
}

genes_up_entrez <- convert_to_entrez(deg_up$gene)
genes_down_entrez <- convert_to_entrez(deg_down$gene)

################################################################################
# 3. GO Enrichment - Upregulated Genes
################################################################################

go_up_bp <- enrichGO(
    gene = genes_up_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP",           # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

go_up_mf <- enrichGO(
    gene = genes_up_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "MF",           # Molecular Function
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
)

# Save results
write.csv(
    as.data.frame(go_up_bp),
    "./results/06_GO/myeloid/GO_upregulated_BP.csv",
    row.names = FALSE
)

write.csv(
    as.data.frame(go_up_mf),
    "./results/06_GO/myeloid/GO_upregulated_MF.csv",
    row.names = FALSE
)

################################################################################
# 4. GO Enrichment - Downregulated Genes
################################################################################

go_down_bp <- enrichGO(
    gene = genes_down_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
)

go_down_mf <- enrichGO(
    gene = genes_down_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
)

write.csv(
    as.data.frame(go_down_bp),
    "./results/06_GO/myeloid/GO_downregulated_BP.csv",
    row.names = FALSE
)

write.csv(
    as.data.frame(go_down_mf),
    "./results/06_GO/myeloid/GO_downregulated_MF.csv",
    row.names = FALSE
)

################################################################################
# 5. Visualizations
################################################################################

# Dot plot - Upregulated BP
if (nrow(as.data.frame(go_up_bp)) > 0) {
    p1 <- dotplot(go_up_bp, showCategory = 20) +
        ggtitle("GO BP - Upregulated in KO")
    
    ggsave(
        "./results/06_GO/myeloid/GO_upregulated_BP_dotplot.png",
        plot = p1, width = 12, height = 8, dpi = 300
    )
}

# Dot plot - Downregulated BP
if (nrow(as.data.frame(go_down_bp)) > 0) {
    p2 <- dotplot(go_down_bp, showCategory = 20) +
        ggtitle("GO BP - Downregulated in KO")
    
    ggsave(
        "./results/06_GO/myeloid/GO_downregulated_BP_dotplot.png",
        plot = p2, width = 12, height = 8, dpi = 300
    )
}

# Bar plot - Top terms
if (nrow(as.data.frame(go_up_bp)) > 0) {
    p3 <- barplot(go_up_bp, showCategory = 15) +
        ggtitle("Top GO Terms - Upregulated")
    
    ggsave(
        "./results/06_GO/myeloid/GO_upregulated_BP_barplot.png",
        plot = p3, width = 10, height = 8, dpi = 300
    )
}

# Clean up
rm(deg_up, deg_down, genes_up_entrez, genes_down_entrez,
   go_up_bp, go_up_mf, go_down_bp, go_down_mf)
gc()

################################################################################
# End of Step 06a
################################################################################
