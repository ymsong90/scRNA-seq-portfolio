################################################################################
# 25-06-19  ▶  CD8 T cell 분석  (porcn.combined.harmony)
#    ① CD8 내부 KO vs WT
#    ② CD8 vs Rest
################################################################################

suppressPackageStartupMessages({
    library(Seurat);      library(dplyr);      library(tibble)
    library(clusterProfiler);  library(org.Mm.eg.db)
    library(enrichplot);  library(ggplot2);    library(stringr)
})

ROOT <- "./figure/250619"
subdir <- function(...){d<-file.path(ROOT, ...); if(!dir.exists(d))dir.create(d,TRUE); d}

## 파라미터 ---------------------------------------------------------------------
PAR <- list(logfc=0.25, minpct=0.10, padj=0.05,
            go_p=0.10, go_q=0.20, ont="BP")

## 공통 GO 실행 함수 ------------------------------------------------------------
run_go <- function(gene_sym, out_tag, out_folder, bg_entrez){
    if(length(gene_sym)<5) return(NULL)
    ent <- bitr(gene_sym,"SYMBOL","ENTREZID",org.Mm.eg.db)$ENTREZID |> unique() |> na.omit()
    if(length(ent)<5) return(NULL)
    ego <- enrichGO(gene=ent, universe=bg_entrez, OrgDb=org.Mm.eg.db,
                    ont=PAR$ont, pvalueCutoff=PAR$go_p, qvalueCutoff=PAR$go_q,
                    readable=TRUE)
    if(nrow(as.data.frame(ego))==0) return(NULL)
    write.csv(as.data.frame(ego), file=file.path(out_folder,paste0(out_tag,".csv")), row.names=FALSE)
    ggsave(file.path(out_folder,paste0(out_tag,"_bar.tiff")),
           barplot(ego,showCategory=10,title=out_tag),
           dpi=600,width=6,height=5,compression="lzw")
    ggsave(file.path(out_folder,paste0(out_tag,"_dot.tiff")),
           dotplot(ego,showCategory=10)+ggtitle(out_tag),
           dpi=600,width=6,height=5,compression="lzw")
    ego
}

## 공통 DEG+GO 워크플로 ---------------------------------------------------------
one_deg_go <- function(seu, ident1, ident2=NULL, prefix, bg_entrez){
    deg <- FindMarkers(seu, ident.1=ident1, ident.2=ident2,
                       logfc.threshold=PAR$logfc, min.pct=PAR$minpct,
                       only.pos=FALSE, verbose=FALSE) |>
        rownames_to_column("gene") |> arrange(p_val_adj)
    # 저장
    write.csv(deg,file=file.path(subdir("DEG_all"),paste0("DEG_",prefix,".csv")),row.names=FALSE)
    up   <- deg %>% filter(avg_log2FC> 0, p_val_adj<PAR$padj)
    down <- deg %>% filter(avg_log2FC< 0, p_val_adj<PAR$padj)
    write.csv(up,  file=file.path(subdir("DEG_up"),  paste0("DEG_Up_",prefix,".csv")),  row.names=FALSE)
    write.csv(down,file=file.path(subdir("DEG_down"),paste0("DEG_Down_",prefix,".csv")),row.names=FALSE)
    # GO
    run_go(up$gene,   paste0(prefix,"_Up"),   subdir("GO_up"),   bg_entrez)
    run_go(down$gene, paste0(prefix,"_Down"), subdir("GO_down"), bg_entrez)
    invisible(list(deg=deg, up=up, down=down))
}

################################################################################
# 분석-1)  CD8 T 내부 KO vs WT
################################################################################
Idents(porcn.combined.harmony) <- "Complete_Labels"
cells_cd8 <- WhichCells(porcn.combined.harmony, idents="CD8 T cell")
cd8_seu   <- subset(porcn.combined.harmony, cells = cells_cd8)

Idents(cd8_seu) <- "orig.ident"    # "porcn_ko" / "porcn_wt"
bg_cd8 <- bitr(rownames(cd8_seu),"SYMBOL","ENTREZID",org.Mm.eg.db)$ENTREZID |> unique() |> na.omit()

one_deg_go(cd8_seu,
           ident1 = "porcn_ko",
           ident2 = "porcn_wt",
           prefix = "CD8_KO_vs_WT",
           bg_entrez = bg_cd8)

################################################################################
# 분석-2)  CD8 T cell  vs  나머지 모든 세포
################################################################################
Idents(porcn.combined.harmony) <- "Complete_Labels"
bg_all <- bitr(rownames(porcn.combined.harmony),"SYMBOL","ENTREZID",org.Mm.eg.db)$ENTREZID |> unique() |> na.omit()

one_deg_go(porcn.combined.harmony,
           ident1 = "CD8 T cell",
           ident2 = NULL,              # rest
           prefix = "CD8_vs_Rest",
           bg_entrez = bg_all)

cat("\n[✓] Finished!  Outputs saved under:", normalizePath(ROOT), "\n")


###############################################################################
# 9) Volcano plots  ────────────────────────────────────────────────────────────
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

volc_dir <- file.path(ROOT, "volcano")
if (!dir.exists(volc_dir)) dir.create(volc_dir, recursive = TRUE)

## ── 재사용할 함수 (컷오프는 PAR$logfc / PAR$padj) ---------------------------
plot_volcano <- function(df, tag,
                         lfc_cut  = PAR$logfc,
                         padj_cut = PAR$padj,
                         label_n  = 20) {          # ← 기본 20
    
    if (nrow(df) == 0) return(NULL)
    
    df <- df %>% mutate(
        log10P   = -log10(p_val_adj + 1e-300),
        sig_flag = case_when(
            p_val_adj < padj_cut & avg_log2FC >  lfc_cut ~ "Up",
            p_val_adj < padj_cut & avg_log2FC < -lfc_cut ~ "Down",
            TRUE ~ "NS"
        ),
        sig_flag = factor(sig_flag, levels = c("Up","Down","NS"))
    )
    
    ## ▲ 여기만 label_n 사용
    top <- df %>% filter(sig_flag != "NS") %>%
        arrange(p_val_adj) %>% group_by(sig_flag) %>%
        slice_head(n = label_n)
    
    ggplot(df, aes(avg_log2FC, log10P, color = sig_flag)) +
        geom_point(size = 4.5, alpha = 0.7) +
        scale_color_manual(breaks=c("Up","Down","NS"),
                           values=c(Up="#f06292",Down="#311b92",NS="grey50")) +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype="dashed") +
        geom_hline(yintercept = -log10(padj_cut),     linetype="dashed") +
        geom_text_repel(data = top, aes(label = gene),
                        size = 4.5, max.overlaps = 50, show.legend = FALSE) +
        labs(title = paste0("Volcano: ", tag),
             x = "log2 Fold Change", y = "-log10(adj.P)") +
        theme_bw() + theme(legend.title = element_blank())
}


## 9-A) CD8 내부 KO vs WT -------------------------------------------------------
raw_cd8_KO_vs_WT <- FindMarkers(cd8_seu,
                                ident.1="porcn_ko", ident.2="porcn_wt",
                                logfc.threshold=0, min.pct=PAR$minpct, verbose=FALSE) |>
    rownames_to_column("gene")
## (A) CD8  KO vs WT  →  기존 컷(0.25)
p1 <- plot_volcano(raw_cd8_KO_vs_WT, "CD8_KO_vs_WT")

ggsave(file.path(volc_dir,"Volcano_CD8_KO_vs_WT.tiff"),
       p1, dpi=600, width=6, height=5, device="tiff", compression="lzw")

## 9-B) CD8 T vs Rest -----------------------------------------------------------
raw_CD8_vs_Rest <- FindMarkers(porcn.combined.harmony,
                               ident.1="CD8 T cell", ident.2=NULL,
                               logfc.threshold=0, min.pct=PAR$minpct, verbose=FALSE) |>
    rownames_to_column("gene")

## (B) CD8 vs Rest ─ logFC 0.5 & 라벨 10개
p2 <- plot_volcano(raw_CD8_vs_Rest,
                   tag      = "CD8_vs_Rest",
                   lfc_cut  = 0.5,
                   label_n  = 10)
ggsave(file.path(volc_dir,"Volcano_CD8_vs_Rest.tiff"),
       p2, dpi=600, width=6, height=5, device="tiff", compression="lzw")

cat("[✓] Volcano plots saved to:", normalizePath(volc_dir), "\n")




################################################################################
# 25-06-19   Mouse CD8 T cells ▶ Dysfunctional / Cytotoxic score + TF correlation
#            (human-signature → mouse ortholog 자동 변환 포함)
################################################################################
#   입력 세션:   cd8_seu  ←  (Seurat v4 object, 이미 CD8 T cell만 추출·정규화 완료)
#   출력 폴더:   ./figure/250619_CD8_DysfCyto_mouse
################################################################################

suppressPackageStartupMessages({
    library(Seurat);    library(dplyr);      library(tibble)
    library(ggplot2);   library(ggrepel);    library(clusterProfiler)
    library(biomaRt)    # ortholog 매핑용
})

# 설치
if (!requireNamespace("homologene", quietly = TRUE))
    install.packages("homologene")

OUT <- "./figure/250619/CD8_DysfCyto_mouse"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1)  Human → Mouse ortholog 매핑 함수
###############################################################################
map_human_to_mouse_local <- function(hsym){
    res <- homologene::homologene(hsym, inTax = 9606, outTax = 10090)
    unique(res$`10090`)
}

###############################################################################
# 2)  Cell 2019 논문 gene-set (human symbol)  →  mouse 기호로 변환
###############################################################################
dysf_human <- c(
    "LAG3","PDCD1","HAVCR2","TIGIT","CXCL13","IFNG","TNFRSF1B","CTLA4","LAYN",
    "ENPP2","BATF","ITM2A","ENTPD1","IL7R","ID3","TOX","TOX2","CXCR6",
    "NR4A1","NR4A2","NR4A3","TNFRSF9","CD38","PTPN23","EOMES","RC3H1","MKI67",
    "TOP2A","HMGB2","CXCL9"
)

cyto_human <- c(
    "FGFBP2","GNLY","CX3CR1","KLRG1","GZMH","PRF1","NKG7","CXCR3","TBX21",
    "GZMB","KLRD1","EOMES","IL7R","IFNG","CCL5","KLF2","DUSP2","GZMA",
    "SELL","GZMK","S1PR1","CD44","CXCR4","ITGAL","ZEB2","STAT1","GBP5",
    "FASLG","RUNX3","IFITM3"
)


library(homologene)
dysf_genes <- map_human_to_mouse_local(dysf_human)
cyto_genes <- map_human_to_mouse_local(cyto_human)

cat("Mapped genes  Dysf:", length(dysf_genes),
    "| Cyto:", length(cyto_genes), "\n")

###############################################################################
# 3)  Module score 계산
###############################################################################
cd8_seu <- AddModuleScore(cd8_seu, features = list(dysf_genes), name = "Dysf")
cd8_seu <- AddModuleScore(cd8_seu, features = list(cyto_genes), name = "Cyto")

###############################################################################
# 4)  2D Scatter 저장
###############################################################################
scatter <- FeatureScatter(cd8_seu, feature1 = "Dysf1", feature2 = "Cyto1",
                          pt.size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    labs(title = "Mouse CD8 T  |  Dysfunctional vs Cytotoxic score",
         x = "Dysfunctional score", y = "Cytotoxic score") +
    theme_bw()
ggsave(file.path(OUT, "Scatter_Dysf_vs_Cyto.tiff"),
       scatter, dpi = 600, width = 6, height = 5, compression = "lzw")

###############################################################################
# 5)  Mouse TF 목록  &  score-correlation
###############################################################################
# Lambert 2018 TF list (human) → mouse로 변환
tf_human <- read.csv("human_TF_lambert2018.csv")$Gene
tf_mouse <- map_human_to_mouse(tf_human)

expr <- GetAssayData(cd8_seu, slot = "data")  # log-normalized
tf_mouse <- intersect(rownames(expr), tf_mouse)

dysf_score <- cd8_seu$Dysf1
cyto_score <- cd8_seu$Cyto1

cor_df <- tibble(
    TF    = tf_mouse,
    r_Dysf = apply(expr[tf_mouse, ], 1, \(x) cor(x, dysf_score, method = "pearson")),
    r_Cyto = apply(expr[tf_mouse, ], 1, \(x) cor(x, cyto_score, method = "pearson"))
)

# 상위 10 TF
top10_dysf <- cor_df %>% arrange(desc(abs(r_Dysf))) %>% slice_head(n = 10)
top10_cyto <- cor_df %>% arrange(desc(abs(r_Cyto))) %>% slice_head(n = 10)

write.csv(top10_dysf, file.path(OUT, "Top10_TF_Dysf.csv"), row.names = FALSE)
write.csv(top10_cyto, file.path(OUT, "Top10_TF_Cyto.csv"), row.names = FALSE)

###############################################################################
# 6)  막대그래프 저장
###############################################################################
plot_bar <- function(df, col, title){
    ggplot(df, aes(x = reorder(TF, !!sym(col)), y = !!sym(col))) +
        geom_col(fill = "#6c8ebf") +
        coord_flip() +
        labs(title = title, x = "", y = "Pearson r") +
        theme_bw()
}

ggsave(file.path(OUT, "Bar_Top10_TF_Dysf.tiff"),
       plot_bar(top10_dysf, "r_Dysf", "Top10 TF vs Dysf score"),
       dpi = 600, width = 5, height = 4, compression = "lzw")

ggsave(file.path(OUT, "Bar_Top10_TF_Cyto.tiff"),
       plot_bar(top10_cyto, "r_Cyto", "Top10 TF vs Cyto score"),
       dpi = 600, width = 5, height = 4, compression = "lzw")

cat("\n[✓] 모든 결과가", normalizePath(OUT), "폴더에 저장되었습니다.\n")



################################################################################
# DotPlot  (Myeloid vs Cancer)  ×  KO / WT  ×  4 chemokines
################################################################################
suppressPackageStartupMessages({
    library(Seurat);  library(dplyr);  library(ggplot2);  library(stringr)
})

out_dir <- "./figure/250619"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

genes_chemokine <- c("Ccl3","Ccl4","Cxcl1","Cxcl2")

# -----------------------------------------------------------------------------#
# 1)  Myeloid 5 클러스터  (WT / KO 별)
# -----------------------------------------------------------------------------#
obj_my <- Myeloid_subset_fil

## 새 그룹 ID  =  "<Cluster>_<WT/KO>"
obj_my$grp <- paste(obj_my$Myeloid_labels, obj_my$orig.ident, sep = "_")
Idents(obj_my) <- "grp"

## DotPlot
p_my <- DotPlot(
    object    = obj_my,
    features  = genes_chemokine,
    cols      = "RdYlBu",
    dot.scale = 10
) + RotatedAxis() +
    theme_bw(base_size = 25) +
    theme(
        axis.text.x  = element_text(angle = 90, vjust = .5, hjust = .5, size = 20),
        axis.text.y  = element_text(size = 20),
        #axis.ticks   = element_blank(),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text  = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.3)
    ) +
    labs(title = "Myeloid clusters  |  KO vs WT", x = NULL, y = NULL)

ggsave(
    filename = file.path(out_dir, "DotPlot_Myeloid_KO_vs_WT.tiff"),
    plot     = p_my,
    dpi      = 600, width = 12, height = 6, compression = "lzw"
)

# -----------------------------------------------------------------------------#
# 2)  Basal / Classical cancer cell  (WT / KO 별)
# -----------------------------------------------------------------------------#
cancer_ids <- c("Basal cancer cell","Classical cancer cell")
obj_ca <- subset(porcn.combined.harmony, idents = cancer_ids)

obj_ca$grp <- paste(obj_ca$Complete_Labels, obj_ca$orig.ident, sep = "_")
Idents(obj_ca) <- "grp"

p_ca <- DotPlot(
    object    = obj_ca,
    features  = genes_chemokine,
    cols      = "RdYlBu",
    dot.scale = 10
) + RotatedAxis() +
    theme_bw(base_size = 25) +
    theme(
        axis.text.x  = element_text(angle = 90, vjust = .5, hjust = .5, size = 20),
        axis.text.y  = element_text(size = 20),
        #axis.ticks   = element_blank(),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text  = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.3)
    ) +
    labs(title = "Cancer cells  |  KO vs WT", x = NULL, y = NULL)

ggsave(
    filename = file.path(out_dir, "DotPlot_Cancer_KO_vs_WT.tiff"),
    plot     = p_ca,
    dpi      = 600, width = 12, height = 6, compression = "lzw"
)

cat("[✓] DotPlots saved to", normalizePath(out_dir), "\n")






