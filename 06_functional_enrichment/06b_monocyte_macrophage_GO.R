#### scRNA-seq Myeloid 세포 DEG 및 GO 분석 (배경유전자 보정 버전) ####

#### 1) 라이브러리 로드 ######################################################
# 필요한 패키지 설치 (최초 1회만)
# BiocManager::install(c("Seurat", "clusterProfiler", "enrichplot", "org.Mm.eg.db", "DOSE"), 
#                     update = TRUE, ask = FALSE)

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)    # 마우스 유전자 매핑용
library(enrichplot)      # 시각화용
library(DOSE)
library(tibble)
library(stringr)
library(ggplot2)

#### 2) 데이터 확인 및 준비 ################################################
# 데이터 상태 확인
cat("=== 데이터 기본 정보 ===\n")
print(dim(Myeloid_subset_fil))
print(table(Myeloid_subset_fil$Myeloid_labels, Myeloid_subset_fil$ID))

# 분석할 클러스터 정의
clusters_to_test <- c(
    "Classical Monocyte",
    "Activated Monocyte", 
    "Hexb+ Macrophage",
    "Selenop+ Macrophage",
    "Trem2+ Macrophage"
)

# 배경 유전자 설정 (매우 중요!)
# scRNA-seq에서는 검출된 유전자만 배경으로 사용해야 함
background_genes <- rownames(Myeloid_subset_fil)
cat("총 검출된 유전자 수:", length(background_genes), "\n")

#### 3) 유전자 ID 매핑 함수 정의 ############################################
# SYMBOL → ENTREZID 변환 함수 (에러 핸들링 포함)
map_symbols_to_entrez <- function(symbols, return_mapping = FALSE) {
    if (length(symbols) == 0) {
        if (return_mapping) {
            return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
        } else {
            return(character(0))
        }
    }
    
    tryCatch({
        mapped <- bitr(
            symbols,
            fromType = "SYMBOL",
            toType = "ENTREZID", 
            OrgDb = org.Mm.eg.db
        )
        
        if (return_mapping) {
            return(mapped)
        } else {
            return(unique(mapped$ENTREZID))
        }
    }, error = function(e) {
        warning("유전자 ID 변환 중 오류: ", e$message)
        if (return_mapping) {
            return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
        } else {
            return(character(0))
        }
    })
}

# 배경 유전자를 ENTREZ ID로 변환
cat("배경 유전자 ID 변환 중...\n")
background_entrez <- map_symbols_to_entrez(background_genes)
cat("변환된 배경 유전자 수:", length(background_entrez), "\n")

#### 4) WT vs KO 차등발현 분석 #############################################
# Seurat ident 설정
Idents(Myeloid_subset_fil) <- "ID"

# 결과 저장용 리스트
deg_list <- list()
deg_summary <- data.frame(
    Cluster = character(),
    Total_DEGs = numeric(),
    Up_in_KO = numeric(), 
    Up_in_WT = numeric(),
    stringsAsFactors = FALSE
)

# 결과 저장 디렉토리 생성
if (!dir.exists("./results")) {
    dir.create("./results", recursive = TRUE)
}

cat("\n=== 차등발현 분석 시작 ===\n")
for (cl in clusters_to_test) {
    cat(">>> 클러스터:", cl, "분석 중...\n")
    
    # 해당 클러스터 셀만 추출
    subset_cl <- subset(Myeloid_subset_fil, subset = Myeloid_labels == cl)
    cell_count <- ncol(subset_cl)
    cat("  - 셀 개수:", cell_count, "\n")
    
    if (cell_count < 10) {
        warning("클러스터 ", cl, "의 셀 수가 너무 적습니다 (n=", cell_count, "). 건너뜁니다.")
        next
    }
    
    # DEG 분석 (더 관대한 기준 적용)
    tryCatch({
        de_markers <- FindMarkers(
            object = subset_cl,
            ident.1 = "KO", 
            ident.2 = "WT",
            logfc.threshold = 0.1,    # 0.25 → 0.1로 완화
            min.pct = 0.05,           # 0.10 → 0.05로 완화  
            test.use = "wilcox",
            verbose = FALSE
        )
        
        # 결과 정리
        de_markers <- de_markers %>%
            rownames_to_column(var = "gene_symbol") %>%
            as_tibble() %>%
            arrange(p_val_adj, desc(abs(avg_log2FC)))
        
        # 유의미한 유전자 필터링
        deg_sig <- de_markers %>% filter(p_val_adj < 0.05)
        
        # 상향/하향 유전자 분류
        up_in_ko <- deg_sig %>% filter(avg_log2FC > 0)
        up_in_wt <- deg_sig %>% filter(avg_log2FC < 0)
        
        cat("  - 총 DEGs:", nrow(deg_sig), 
            " (KO 상향:", nrow(up_in_ko), 
            ", WT 상향:", nrow(up_in_wt), ")\n")
        
        # 결과 저장
        deg_list[[cl]] <- deg_sig
        
        # 요약 정보 추가
        deg_summary <- rbind(deg_summary, data.frame(
            Cluster = cl,
            Total_DEGs = nrow(deg_sig),
            Up_in_KO = nrow(up_in_ko),
            Up_in_WT = nrow(up_in_wt)
        ))
        
        # CSV 저장
        write.csv(
            deg_sig,
            file = paste0("./results/DEG_", 
                          str_replace_all(cl, "\\+", "plus"), 
                          "_KO_vs_WT.csv"),
            row.names = FALSE
        )
        
    }, error = function(e) {
        warning("클러스터 ", cl, "에서 DEG 분석 실패: ", e$message)
    })
}

# DEG 요약 출력 및 저장
cat("\n=== DEG 분석 요약 ===\n")
print(deg_summary)
write.csv(deg_summary, "./results/DEG_summary.csv", row.names = FALSE)

#### 5) GO 분석 함수 정의 ##################################################
perform_go_analysis <- function(gene_list, universe_genes, analysis_name, 
                                ont_type = "BP", pval_cutoff = 0.1, qval_cutoff = 0.2) {
    
    if (length(gene_list) == 0) {
        warning(analysis_name, ": 분석할 유전자가 없습니다.")
        return(NULL)
    }
    
    # 유전자 ID 변환
    entrez_genes <- map_symbols_to_entrez(gene_list)
    entrez_universe <- map_symbols_to_entrez(universe_genes)
    
    cat("  - 입력 유전자:", length(gene_list), 
        "→ 변환된 ENTREZ ID:", length(entrez_genes), "\n")
    
    if (length(entrez_genes) < 5) {
        warning(analysis_name, ": 변환된 유전자 수가 너무 적습니다 (n=", 
                length(entrez_genes), "). 건너뜁니다.")
        return(NULL)
    }
    
    # GO 분석 실행
    tryCatch({
        ego <- enrichGO(
            gene = entrez_genes,
            universe = entrez_universe,  # 배경 유전자 설정 (핵심!)
            OrgDb = org.Mm.eg.db,
            keyType = "ENTREZID",
            ont = ont_type,
            pAdjustMethod = "BH", 
            pvalueCutoff = pval_cutoff,
            qvalueCutoff = qval_cutoff,
            readable = TRUE
        )
        
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
            cat("  - GO 분석 성공:", nrow(as.data.frame(ego)), "개 GO term 발견\n")
            return(ego)
        } else {
            warning(analysis_name, ": 유의미한 GO term이 없습니다.")
            return(NULL)
        }
        
    }, error = function(e) {
        warning(analysis_name, ": GO 분석 실패 - ", e$message)
        return(NULL)
    })
}

#### 6) GO 분석 실행 #######################################################
go_results_ko_up <- list()  # KO에서 상향된 유전자
go_results_wt_up <- list()  # WT에서 상향된 유전자

cat("\n=== GO 분석 시작 ===\n")
for (cl in names(deg_list)) {
    cat(">>> 클러스터:", cl, "GO 분석 중...\n")
    
    deg_data <- deg_list[[cl]]
    
    if (nrow(deg_data) == 0) {
        cat("  - DEG가 없어 건너뜁니다.\n")
        next
    }
    
    # KO에서 상향된 유전자 분석
    genes_up_ko <- deg_data %>% filter(avg_log2FC > 0) %>% pull(gene_symbol)
    if (length(genes_up_ko) > 0) {
        cat("  [KO 상향 유전자 GO 분석]\n")
        go_ko <- perform_go_analysis(
            gene_list = genes_up_ko,
            universe_genes = background_genes,
            analysis_name = paste(cl, "- KO up")
        )
        go_results_ko_up[[cl]] <- go_ko
        
        # 결과 저장
        if (!is.null(go_ko)) {
            write.csv(
                as.data.frame(go_ko),
                file = paste0("./results/GO_BP_", 
                              str_replace_all(cl, "\\+", "plus"), 
                              "_KO_up.csv"),
                row.names = FALSE
            )
        }
    }
    
    # WT에서 상향된 유전자 분석  
    genes_up_wt <- deg_data %>% filter(avg_log2FC < 0) %>% pull(gene_symbol)
    if (length(genes_up_wt) > 0) {
        cat("  [WT 상향 유전자 GO 분석]\n")
        go_wt <- perform_go_analysis(
            gene_list = genes_up_wt,
            universe_genes = background_genes,
            analysis_name = paste(cl, "- WT up")
        )
        go_results_wt_up[[cl]] <- go_wt
        
        # 결과 저장
        if (!is.null(go_wt)) {
            write.csv(
                as.data.frame(go_wt),
                file = paste0("./results/GO_BP_", 
                              str_replace_all(cl, "\\+", "plus"), 
                              "_WT_up.csv"),
                row.names = FALSE
            )
        }
    }
    
    cat("\n")
}

#### 7) 결과 시각화 ########################################################
# 시각화 함수
create_go_plots <- function(go_result, title_prefix, cluster_name) {
    if (is.null(go_result) || nrow(as.data.frame(go_result)) == 0) {
        return(NULL)
    }
    
    plots <- list()
    
    # Bar plot
    tryCatch({
        p1 <- barplot(go_result, showCategory = 10, 
                      title = paste(title_prefix, cluster_name, "GO:BP (Top 10)"))
        plots$barplot <- p1
    }, error = function(e) {
        warning("Bar plot 생성 실패: ", e$message)
    })
    
    # Dot plot
    tryCatch({
        p2 <- dotplot(go_result, showCategory = 10) + 
            ggtitle(paste(title_prefix, cluster_name, "GO:BP")) +
            theme(axis.text.y = element_text(size = 8))
        plots$dotplot <- p2
    }, error = function(e) {
        warning("Dot plot 생성 실패: ", e$message)
    })
    
    return(plots)
}

# 시각화 실행 및 개별 TIFF 파일로 저장
cat("=== 시각화 생성 중 ===\n")

# 저장 디렉토리 확인 및 생성
if (!dir.exists("./results/GO_plots")) {
    dir.create("./results/GO_plots", recursive = TRUE)
}

for (cl in names(go_results_ko_up)) {
    # 클러스터 이름에서 파일명에 사용할 수 없는 문자 제거
    clean_cluster_name <- str_replace_all(cl, "\\+", "plus")
    clean_cluster_name <- str_replace_all(clean_cluster_name, " ", "_")
    
    # KO 상향 유전자 플롯
    if (!is.null(go_results_ko_up[[cl]])) {
        cat("KO 상향 -", cl, "플롯 생성 중...\n")
        plots_ko <- create_go_plots(go_results_ko_up[[cl]], "KO Up-regulated", cl)
        
        # Bar plot 저장
        if (!is.null(plots_ko$barplot)) {
            ggsave(
                filename = paste0("./results/GO_plots/250602_GO_barplot_", clean_cluster_name, "_KO_up.tiff"),
                plot     = plots_ko$barplot,
                device   = "tiff",
                dpi      = 600,
                width    = 12,
                height   = 8,
                units    = "in",
                compression = "lzw"
            )
        }
        
        # Dot plot 저장
        if (!is.null(plots_ko$dotplot)) {
            ggsave(
                filename = paste0("./results/GO_plots/250602_GO_dotplot_", clean_cluster_name, "_KO_up.tiff"),
                plot     = plots_ko$dotplot,
                device   = "tiff",
                dpi      = 600,
                width    = 12,
                height   = 8,
                units    = "in",
                compression = "lzw"
            )
        }
    }
    
    # WT 상향 유전자 플롯
    if (!is.null(go_results_wt_up[[cl]])) {
        cat("WT 상향 -", cl, "플롯 생성 중...\n")
        plots_wt <- create_go_plots(go_results_wt_up[[cl]], "WT Up-regulated", cl)
        
        # Bar plot 저장
        if (!is.null(plots_wt$barplot)) {
            ggsave(
                filename = paste0("./results/GO_plots/250602_GO_barplot_", clean_cluster_name, "_WT_up.tiff"),
                plot     = plots_wt$barplot,
                device   = "tiff",
                dpi      = 600,
                width    = 12,
                height   = 8,
                units    = "in",
                compression = "lzw"
            )
        }
        
        # Dot plot 저장
        if (!is.null(plots_wt$dotplot)) {
            ggsave(
                filename = paste0("./results/GO_plots/250602_GO_dotplot_", clean_cluster_name, "_WT_up.tiff"),
                plot     = plots_wt$dotplot,
                device   = "tiff",
                dpi      = 600,
                width    = 12,
                height   = 8,
                units    = "in",
                compression = "lzw"
            )
        }
    }
}

cat("모든 GO 플롯이 ./results/GO_plots/ 폴더에 TIFF 형식으로 저장되었습니다.\n")

#### 8) 분석 결과 요약 ####################################################
cat("\n=== 최종 분석 결과 요약 ===\n")

# GO 분석 요약
go_summary <- data.frame(
    Cluster = character(),
    KO_up_GO_terms = numeric(),
    WT_up_GO_terms = numeric(),
    stringsAsFactors = FALSE
)

for (cl in clusters_to_test) {
    ko_terms <- if (!is.null(go_results_ko_up[[cl]])) nrow(as.data.frame(go_results_ko_up[[cl]])) else 0
    wt_terms <- if (!is.null(go_results_wt_up[[cl]])) nrow(as.data.frame(go_results_wt_up[[cl]])) else 0
    
    go_summary <- rbind(go_summary, data.frame(
        Cluster = cl,
        KO_up_GO_terms = ko_terms,
        WT_up_GO_terms = wt_terms
    ))
}

print("DEG 요약:")
print(deg_summary)
cat("\nGO 분석 요약:\n")
print(go_summary)

# 최종 요약 저장
write.csv(go_summary, "./results/250602_GO_summary.csv", row.names = FALSE)

cat("\n분석 완료! 결과는 ./results/ 폴더에 저장되었습니다.\n")
cat("- DEG 결과: DEG_[클러스터명]_KO_vs_WT.csv\n")
cat("- GO 결과: GO_BP_[클러스터명]_[KO/WT]_up.csv\n")
cat("- 시각화: GO_analysis_plots.pdf\n")
cat("- 요약: DEG_summary.csv, GO_summary.csv\n")
