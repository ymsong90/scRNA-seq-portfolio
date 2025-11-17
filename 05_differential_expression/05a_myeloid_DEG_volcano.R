# -----------------------------------------------------------------------------
# 1) 클러스터별 (Myeloid_Lables) WT vs KO DEG 분석
# -----------------------------------------------------------------------------
Idents(Myeloid_subset_fil) <- "Myeloid_labels"

# “클러스터+조건”이라는 새로운 메타데이터 생성
Myeloid_subset_fil$cluster.condition <- paste0(
    Idents(Myeloid_subset_fil), "_", Myeloid_subset_fil$orig.ident
)

Idents(Myeloid_subset_fil) <- "cluster.condition"
cluster_ids <- unique(stringr::str_split_fixed(
    levels(Idents(Myeloid_subset_fil)), "_", 2
)[,1])

# 각 클러스터별 FindMarkers 수행, 그리고 rownames(=gene) → gene 컬럼으로
cluster_degs <- lapply(cluster_ids, function(cl){
    ident1 <- paste0(cl, "_porcn_ko")
    ident2 <- paste0(cl, "_porcn_wt")
    if (!(ident1 %in% levels(Idents(Myeloid_subset_fil))) ||
        !(ident2 %in% levels(Idents(Myeloid_subset_fil))))
        return(NULL)
    fm <- FindMarkers(
        Myeloid_subset_fil,
        ident.1         = ident1,
        ident.2         = ident2,
        min.pct         = 0.1,
        logfc.threshold = 0.25,
        verbose         = FALSE
    )
    # gene 이름을 컬럼으로 끌어오기
    fm <- rownames_to_column(fm, var = "gene")
    fm$cluster <- cl
    fm
})

names(cluster_degs) <- cluster_ids

# 클러스터 3 예시 확인
head(cluster_degs[["Hexb+ Macrophage"]], 10)

# 모든 클러스터 결과 합치기 및 저장
all_cluster_degs <- bind_rows(cluster_degs, .id = "cluster")
write.csv(all_cluster_degs,
          file      = "./figure/250618/250618_MonoMac_perCluster_DEG_ko_vs_wt.csv",
          row.names = FALSE)




################################################################################
#  Myeloid scRNA-seq ▶ DEG (KO vs WT) ▶ GO (Up & Down) ▶ Plot / Export
#  Last updated: 2025-06-18
################################################################################

## 0) 환경설정 ------------------------------------------------------------------
suppressPackageStartupMessages({
    library(Seurat);            library(dplyr);       library(tibble)
    library(clusterProfiler);   library(org.Mm.eg.db)
    library(enrichplot);        library(ggplot2);     library(stringr)
})

# ── (선택) 병렬 처리 ----------------------------------------------------------
if (requireNamespace("future.apply", quietly = TRUE)) {
    library(future.apply)
    plan(multisession, workers = max(1, parallel::detectCores() - 1))
}

## 1) 전역 파라미터 & 유틸 ------------------------------------------------------
PARAMS <- list(
    outdir      = "./figure/250618",
    deg         = list(logfc = 0.25, minpct = 0.10, padj = 0.05),
    go          = list(pval = 0.10, qval = 0.20, ont = "BP", showN = 10)
)

create_dir <- function(p) if (!dir.exists(p)) dir.create(p, TRUE)
dirs <- c(PARAMS$outdir,
          file.path(PARAMS$outdir, "DEG"),
          file.path(PARAMS$outdir, "GO"),
          file.path(PARAMS$outdir, "plots"))
invisible(lapply(dirs, create_dir))

sanitize <- function(x) str_replace_all(x, "[^A-Za-z0-9]", "_")

## 2) 배경 유전자 ---------------------------------------------------------------
background_entrez_ids <- bitr(
    rownames(Myeloid_subset_fil),
    fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db
)$ENTREZID |> unique()
message("[Info] Background genes mapped: ", length(background_entrez_ids))

## 3) 클러스터별 DEG (raw + filtered) ------------------------------------------
Idents(Myeloid_subset_fil) <- "Myeloid_labels"
Myeloid_subset_fil$cluster.condition <- paste0(
    Idents(Myeloid_subset_fil), "_", Myeloid_subset_fil$orig.ident)
Idents(Myeloid_subset_fil) <- "cluster.condition"
cluster_ids <- unique(str_split_fixed(levels(Idents(Myeloid_subset_fil)), "_", 2)[,1])

# 3-A) RAW -- 모든 유전자 (logFC 0, p-adj 필터 X) ------------------------------
cluster_raw <- future_lapply(cluster_ids, function(cl){
    id1 <- paste0(cl, "_porcn_ko")
    id2 <- paste0(cl, "_porcn_wt")
    if (!(id1 %in% levels(Idents(Myeloid_subset_fil))) ||
        !(id2 %in% levels(Idents(Myeloid_subset_fil)))) return(NULL)
    
    FindMarkers(Myeloid_subset_fil,
                ident.1 = id1, ident.2 = id2,
                logfc.threshold = 0,               # ← 필터 해제
                min.pct = PARAMS$deg$minpct,
                verbose = FALSE) |>
        rownames_to_column("gene") |> mutate(cluster = cl)
}, future.seed = TRUE)
names(cluster_raw) <- cluster_ids   # 리스트 이름 지정

# 3-B) FILTERED -- 기존 컷오프 적용 -------------------------------------------
cluster_degs <- lapply(cluster_raw, function(df){
    if (is.null(df)) return(NULL)
    df %>% filter(abs(avg_log2FC) > PARAMS$deg$logfc,
                  p_val_adj       < PARAMS$deg$padj)
})
write.csv(bind_rows(cluster_degs, .id = "cluster"),
          file = file.path(PARAMS$outdir, "DEG", "Myeloid_DEG_ko_vs_wt.csv"),
          row.names = FALSE)

## 4) GO 함수 -------------------------------------------------------------------
perform_go <- function(genes){
    if (length(genes) < 5) return(NULL)
    entrez <- bitr(genes, "SYMBOL", "ENTREZID", org.Mm.eg.db)$ENTREZID |> unique()
    if (length(entrez) < 5) return(NULL)
    enrichGO(gene=entrez, universe=background_entrez_ids,
             OrgDb=org.Mm.eg.db, ont=PARAMS$go$ont,
             pAdjustMethod="BH",
             pvalueCutoff=PARAMS$go$pval, qvalueCutoff=PARAMS$go$qval,
             readable=TRUE)
}

save_barplot <- function(go_obj, tag){
    if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) return()
    ggsave(
        file.path(PARAMS$outdir, "plots", paste0(tag, ".tiff")),
        barplot(go_obj, showCategory = PARAMS$go$showN, title = tag),
        dpi = 600, width = 10, height = 6, units = "in",
        device = "tiff", compression = "lzw"
    )
}

## 5) GO (KO_up & KO_down) ------------------------------------------------------
go_summary <- tibble(Cluster=character(), KO_up_GO=integer(), KO_down_GO=integer())

for (cl in names(cluster_degs)) {
    dd <- cluster_degs[[cl]]; if (is.null(dd) || nrow(dd)==0) next
    
    # ── KO ↑ (avg_log2FC > 0) ----------------------------------------------------
    ko_up <- dd %>% filter(avg_log2FC > 0) %>% pull(gene)
    go_up <- perform_go(ko_up)
    if (!is.null(go_up)){
        tag <- paste0("GO_", sanitize(cl), "_KO_up")
        write.csv(as.data.frame(go_up),
                  file = file.path(PARAMS$outdir, "GO", paste0(tag, ".csv")),
                  row.names = FALSE)
        save_barplot(go_up, tag)
    }
    
    # ── KO ↓ (avg_log2FC < 0)  (= WT ↑) ----------------------------------------
    ko_down <- dd %>% filter(avg_log2FC < 0) %>% pull(gene)
    go_down <- perform_go(ko_down)
    if (!is.null(go_down)){
        tag <- paste0("GO_", sanitize(cl), "_KO_down")
        write.csv(as.data.frame(go_down),
                  file = file.path(PARAMS$outdir, "GO", paste0(tag, ".csv")),
                  row.names = FALSE)
        save_barplot(go_down, tag)
    }
    
    go_summary <- add_row(go_summary,
                          Cluster      = cl,
                          KO_up_GO     = ifelse(is.null(go_up),   0, nrow(go_up)),
                          KO_down_GO   = ifelse(is.null(go_down), 0, nrow(go_down)))
}

write.csv(go_summary, file.path(PARAMS$outdir, "GO_summary.csv"), row.names = FALSE)
message("\n[✓] Pipeline finished!  All results are under: ", normalizePath(PARAMS$outdir))


## 6) Volcano plot  ────────────────────────────────────────────────────────────
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

volc_dir <- file.path(PARAMS$outdir, "plots", "volcano")
create_dir(volc_dir)

# ── 시각화 함수 ---------------------------------------------------------------
plot_volcano <- function(df, cl){
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    volc_cut <- list(p_adj = PARAMS$deg$padj, logfc = PARAMS$deg$logfc)
    
    df <- df %>%
        mutate(
            log10P   = -log10(p_val_adj + 1e-300),
            sig_flag = case_when(
                p_val_adj < volc_cut$p_adj & avg_log2FC >  volc_cut$logfc ~ "Up",
                p_val_adj < volc_cut$p_adj & avg_log2FC < -volc_cut$logfc ~ "Down",
                TRUE ~ "NS"
            ),
            sig_flag = factor(sig_flag, levels = c("Up","Down","NS"))   # ② 순서 고정
        )
    
    top_labs <- df %>% filter(sig_flag != "NS") %>% 
        arrange(p_val_adj) %>% group_by(sig_flag) %>% slice_head(n = 20)
    
    ggplot(df, aes(avg_log2FC, log10P, color = sig_flag)) +
        geom_point(size = 2, alpha = 0.7) +
        scale_color_manual(
            breaks = c("Up","Down","NS"),                               # ②
            values = c(Up="#f06292", Down="#311b92", NS="grey50")
        ) +
        geom_vline(xintercept = c(-volc_cut$logfc, volc_cut$logfc), linetype="dashed") +
        geom_hline(yintercept = -log10(volc_cut$p_adj), linetype="dashed") +
        geom_text_repel(                                           # ① legend 끔
            data = top_labs, aes(label = gene),
            size = 3, max.overlaps = 50, show.legend = FALSE
        ) +
        labs(title = paste0("Volcano: ", cl, " (KO vs WT)"),
             x = "log2 Fold Change", y = "-log10(adj.P)") +
        theme_bw() + theme(legend.title = element_blank())
}

# ── 루프: 반드시 cluster_raw 사용 -------------------------------------------
for (cl in names(cluster_raw)) {
    raw_df <- cluster_raw[[cl]]
    p <- plot_volcano(raw_df, cl)
    if (!is.null(p)) {
        ggsave(
            filename = file.path(volc_dir, paste0("Volcano_", sanitize(cl), ".tiff")),
            plot     = p,
            dpi      = 600, width = 6, height = 5, units = "in",
            device   = "tiff", compression = "lzw"
        )
    }
}



###############################################################################
# 25-06-18. DEG_one_vs_rest, GO_one_vs_rest 분석 (Myeloid_subset_fil)
# (1) 필수 패키지 로드  (동일)
###############################################################################
suppressPackageStartupMessages({
    library(Seurat);   library(dplyr);   library(tibble)
    library(clusterProfiler);  library(org.Mm.eg.db)
    library(enrichplot);  library(ggplot2);  library(stringr)
})

###############################################################################
# (2) 출력 폴더 (Up / Down 별 하위폴더 추가)   ◆ NEW
###############################################################################
out_deg <- "./figure/250618/DEG_one_vs_rest"
out_go  <- "./figure/250618/GO_one_vs_rest"
dir.create(file.path(out_deg, "Up"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_deg, "Down"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_go,  "Up"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_go,  "Down"), recursive = TRUE, showWarnings = FALSE)

###############################################################################
# (3) 클러스터 레벨 정의 (동일)
###############################################################################
cl_levels <- levels(factor(Myeloid_subset_fil$Myeloid_labels))

###############################################################################
# (4) 배경 ENTREZ ID 셋 (동일)
###############################################################################
bg_symbols <- rownames(Myeloid_subset_fil)
bg_entrez  <- bitr(bg_symbols, "SYMBOL", "ENTREZID", org.Mm.eg.db)$ENTREZID |>
    unique() |>
    na.omit()

###############################################################################
# (5) One-vs-Rest  ▶  Up / Down  DEG + GO
###############################################################################
Idents(Myeloid_subset_fil) <- "Myeloid_labels" 

deg_summary <- list()

for (cl in cl_levels) {
    
    message("── ", cl, "  (one-vs-rest)")
    
    ## 5-1  DEG  ---------------------------------------------------------------
    deg <- FindMarkers(
        object          = Myeloid_subset_fil,
        ident.1         = cl,
        ident.2         = NULL,          # 나머지 전체
        only.pos        = FALSE,
        logfc.threshold = 0.25,
        min.pct         = 0.10
    ) |>
        rownames_to_column("gene") |>
        arrange(p_val_adj)
    
    ## 5-1-a  Up / Down 분리   ◆ NEW
    deg_up   <- deg %>% filter(avg_log2FC >  0, p_val_adj < 0.05)
    deg_down <- deg %>% filter(avg_log2FC <  0, p_val_adj < 0.05)
    
    write.csv(deg_up,
              file = file.path(out_deg, "Up",
                               paste0("DEG_Up_", str_replace_all(cl, "[ +]", "_"), "_vs_rest.csv")),
              row.names = FALSE)
    write.csv(deg_down,
              file = file.path(out_deg, "Down",
                               paste0("DEG_Down_", str_replace_all(cl, "[ +]", "_"), "_vs_rest.csv")),
              row.names = FALSE)
    
    deg_summary[[cl]] <- c(Up = nrow(deg_up), Down = nrow(deg_down))
    
    ## 5-2  공통 GO 실행 함수  --------------------------------------------------
    run_go <- function(sym_vec, direction){
        if (length(sym_vec) < 5) return(NULL)     # 최소 5개
        ent <- bitr(sym_vec, "SYMBOL", "ENTREZID", org.Mm.eg.db)$ENTREZID |>
            unique() |>
            na.omit()
        if (length(ent) < 5) return(NULL)
        
        ego <- enrichGO(
            gene     = ent,
            universe = bg_entrez,
            OrgDb    = org.Mm.eg.db,
            ont      = "BP",
            pvalueCutoff = 0.10,
            qvalueCutoff = 0.20,
            readable = TRUE
        )
        if (nrow(as.data.frame(ego)) == 0) return(NULL)
        
        tag <- paste0("GO_BP_", str_replace_all(cl, "[ +]", "_"), "_", direction)
        write.csv(as.data.frame(ego),
                  file = file.path(out_go, direction, paste0(tag, ".csv")),
                  row.names = FALSE)
        
        # Barplot & Dotplot (top 10)
        barp <- barplot(ego, showCategory = 10,
                        title = paste0(cl, " (", direction, ")")) +
            theme(plot.title = element_text(size = 14, face = "bold"))
        dotp <- dotplot(ego, showCategory = 10, font.size = 10) +
            ggtitle(paste0(cl, " (", direction, ")"))
        
        ggsave(file.path(out_go, direction, paste0(tag, "_bar.tiff")),
               barp, dpi=600, width=6, height=5, compression="lzw")
        ggsave(file.path(out_go, direction, paste0(tag, "_dot.tiff")),
               dotp, dpi=600, width=6, height=5, compression="lzw")
    }
    
    ## 5-3  GO(Up) & GO(Down)  ◆ NEW
    run_go(deg_up$gene,   "Up")
    run_go(deg_down$gene, "Down")
}

###############################################################################
# (6) 요약 테이블 저장  (Up / Down 열 추가)  ◆ NEW
###############################################################################
deg_summary_df <- tibble(
    Cluster = names(deg_summary),
    Up_DEG  = sapply(deg_summary, `[[`, "Up"),
    Down_DEG= sapply(deg_summary, `[[`, "Down")
)
write.csv(deg_summary_df,
          "./results/DEG_one_vs_rest_summary.csv",
          row.names = FALSE)


###############################################################################
# 25-06-18. DEG_one_vs_rest, GO_one_vs_rest 분석 (porcn.combined.harmony 클러스터 나누지 않은 전체 세포에 대해서서)
###############################################################################

suppressPackageStartupMessages({
    library(Seurat);  library(dplyr);  library(tibble)
    library(clusterProfiler);  library(org.Mm.eg.db)
    library(enrichplot);  library(ggplot2);  library(stringr)
})

###############################################################################
# 1) 출력 폴더
###############################################################################
out_dir <- "./figure/250618/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2) Idents를 KO / WT 로 설정
###############################################################################
Idents(porcn.combined.harmony) <- "orig.ident"     # 예: "porcn_ko" / "porcn_wt"

###############################################################################
# 3) 배경 유전자
###############################################################################
bg_entrez <- bitr(rownames(porcn.combined.harmony), "SYMBOL", "ENTREZID",
                  org.Mm.eg.db)$ENTREZID |> unique() |> na.omit()

###############################################################################
# 4) DEG (전체 세포)
###############################################################################
deg <- FindMarkers(
    porcn.combined.harmony,
    ident.1         = "porcn_ko",
    ident.2         = "porcn_wt",
    only.pos        = FALSE,
    logfc.threshold = 0.25,
    min.pct         = 0.10
) |>
    rownames_to_column("gene") |>
    arrange(p_val_adj)

write.csv(deg, file = file.path(out_dir, "DEG_total_KO_vs_WT.csv"), row.names = FALSE)

###############################################################################
# 5) Up / Down 분리
###############################################################################
up_genes   <- deg %>% filter(avg_log2FC >  0, p_val_adj < 0.05) %>% pull(gene)
down_genes <- deg %>% filter(avg_log2FC <  0, p_val_adj < 0.05) %>% pull(gene)

###############################################################################
# 6) GO 함수
###############################################################################
run_go <- function(sym_vec, tag){
    if (length(sym_vec) < 5) return(NULL)
    ent <- bitr(sym_vec, "SYMBOL", "ENTREZID", org.Mm.eg.db)$ENTREZID |>
        unique() |> na.omit()
    if (length(ent) < 5) return(NULL)
    ego <- enrichGO(
        gene      = ent,
        universe  = bg_entrez,
        OrgDb     = org.Mm.eg.db,
        ont       = "BP",
        pvalueCutoff = 0.10,
        qvalueCutoff = 0.20,
        readable  = TRUE
    )
    if (nrow(as.data.frame(ego)) == 0) return(NULL)
    
    write.csv(as.data.frame(ego),
              file = file.path(out_dir, paste0("GO_BP_total_", tag, ".csv")),
              row.names = FALSE)
    
    ggsave(file.path(out_dir, paste0("GO_bar_total_", tag, ".tiff")),
           barplot(ego, showCategory = 10, title = tag),
           dpi = 600, width = 6, height = 5, compression = "lzw")
    invisible(ego)
}

go_up   <- run_go(up_genes,   "KO_up")
go_down <- run_go(down_genes, "KO_down")

###############################################################################
# 7) 간단 요약
###############################################################################
cat("Up-DEG :", length(up_genes), "\nDown-DEG :", length(down_genes), "\n")




















