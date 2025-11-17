###############################################################################
# (1) 필수 패키지 로드  ─ tibble 포함  +  경고 억제 옵션
###############################################################################
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tibble)   # ← rownames_to_column(), column_to_rownames() 제공
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(stringr)
})

###############################################################################
# (2) 출력 폴더 준비  ─ 없으면 생성
###############################################################################
out_deg <- "./results/DEG_one_vs_rest"
out_go  <- "./results/GO_one_vs_rest"
dir.create(out_deg, showWarnings = FALSE, recursive = TRUE)
dir.create(out_go,  showWarnings = FALSE, recursive = TRUE)

###############################################################################
# (3) 클러스터 레벨 정의
###############################################################################
cl_levels <- levels(factor(Myeloid_subset_fil$Myeloid_labels))
# 또는 원하는 순서가 있다면 ← 수동 지정
# cl_levels <- c("Trem2+ Macrophage","Selenop+ Macrophage",
#                "Hexb+ Macrophage","Activated Monocyte",
#                "Classical Monocyte")

###############################################################################
# (4) 배경 ENTREZ ID 셋 만들기 (NA 제거·중복 제거)
###############################################################################
bg_symbols <- rownames(Myeloid_subset_fil)
bg_entrez  <- bitr(bg_symbols,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Mm.eg.db) |>
    dplyr::pull(ENTREZID) |>
    unique() |>
    na.omit()

###############################################################################
# (5) One-vs-Rest DEG + GO 루프
###############################################################################
deg_summary1 <- list()

for (cl in cl_levels) {
    
    message("── ", cl, "  (one-vs-rest)")
    
    ## 5-1 DEG ------------------------------------------------------------------
    deg <- FindMarkers(
        object          = Myeloid_subset_fil,
        ident.1         = cl,
        ident.2         = NULL,              # 나머지 전체
        only.pos        = FALSE,
        logfc.threshold = 0.25,
        min.pct         = 0.10
    ) |>
        rownames_to_column(var = "gene") |>
        arrange(p_val_adj)
    
    write.csv(
        deg,
        file = file.path(out_deg,
                         paste0("DEGs_", str_replace_all(cl, "[ +]", "_"),
                                "_vs_rest.csv")),
        row.names = FALSE
    )
    
    deg_summary1[[cl]] <- nrow(filter(deg, p_val_adj < 0.05))
    
    ## 5-2 GO (상향 gene만) ------------------------------------------------------
    up_sym <- deg |>
        filter(avg_log2FC > 0, p_val_adj < 0.05) |>
        pull(gene)
    
    if (length(up_sym) >= 5) {
        
        up_ent <- bitr(up_sym,
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Mm.eg.db) |>
            pull(ENTREZID) |>
            unique() |>
            na.omit()
        
        ego <- enrichGO(
            gene      = up_ent,
            universe  = bg_entrez,
            OrgDb     = org.Mm.eg.db,
            ont       = "BP",
            pvalueCutoff = 0.10,
            qvalueCutoff = 0.20,
            readable  = TRUE
        )
        
        if (nrow(as.data.frame(ego)) > 0) {
            
            ## CSV 저장
            write.csv(
                as.data.frame(ego),
                file = file.path(out_go,
                                 paste0("GO_BP_",
                                        str_replace_all(cl, "[ +]", "_"),
                                        "_up.csv")),
                row.names = FALSE
            )
            
            ## 시각화 (top 10)
            barp <- barplot(ego, showCategory = 10,
                            title = paste(cl, "(up-genes)")) +
                theme(plot.title = element_text(size = 14, face = "bold"))
            
            dotp <- dotplot(ego, showCategory = 10, font.size = 10) +
                ggtitle(paste(cl, "(up-genes)"))
            
            ggsave(
                filename   = file.path(out_go,
                                       paste0("GO_bar_",
                                              str_replace_all(cl, "[ +]", "_"),
                                              ".tiff")),
                plot       = barp,
                dpi        = 600,
                width      = 6,
                height     = 5,
                compression = "lzw"
            )
            
            ggsave(
                filename   = file.path(out_go,
                                       paste0("GO_dot_",
                                              str_replace_all(cl, "[ +]", "_"),
                                              ".tiff")),
                plot       = dotp,
                dpi        = 600,
                width      = 6,
                height     = 5,
                compression = "lzw"
            )
        }
    }
}

###############################################################################
# (6) 요약 테이블 저장
###############################################################################
deg_summary1_df <- tibble(
    Cluster        = names(deg_summary1),
    DEG_padj_lt_0.05 = unlist(deg_summary1)
)
write.csv(deg_summary1_df,
          "./results/DEG_one_vs_rest_summary.csv",
          row.names = FALSE)





