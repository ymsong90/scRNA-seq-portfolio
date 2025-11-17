## "Seurat v5 호환 Wnt 신호전달 경로 유전자 발현 분석 및 골수세포 아형별 시각화"


# ──────────────────────────────────────────────────────────
# 0) 패키지
# ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr); library(tidyr)
    library(ggplot2); library(ggthemes); library(forcats)
    library(Matrix)     # sparse 연산
})

# ──────────────────────────────────────────────────────────
# 1) 유전자 목록 (Wnt 리간드 추가)
# ──────────────────────────────────────────────────────────
# Wnt 리간드
wnt_features <- c("Wnt1", "Wnt2", "Wnt2b", "Wnt3", "Wnt4", "Wnt5a", "Wnt5b", "Wnt6",
                  "Wnt7a", "Wnt7b", "Wnt8a", "Wnt8b", "Wnt9a", "Wnt9b", "Wnt10a", "Wnt10b",
                  "Wnt11", "Wnt16")

# 리셉터
rct_features <- c('Lrp1','Lrp5','Lrp6',
                  'Fzd1','Fzd2','Fzd3','Fzd4','Fzd5',
                  'Fzd6','Ror1','Ror2','Ryk')

# 전사인자
tf_features  <- c('Tcf3','Tcf4','Tcf7','Lef1')

# 전체 유전자 세트
gene_set     <- unique(c(wnt_features, rct_features, tf_features))

# ──────────────────────────────────────────────────────────
# 2) WT 세포 서브셋
# ──────────────────────────────────────────────────────────
my_wt <- subset(Myeloid_subset_fil, subset = ID == "WT")
DefaultAssay(my_wt) <- "RNA"

# ──────────────────────────────────────────────────────────
# 3) Seurat v5 호환성 확인 및 설정
# ──────────────────────────────────────────────────────────
# 현재 레이어 상태 확인
cat("현재 레이어:", paste(Layers(my_wt[["RNA"]]), collapse = ", "), "\n")

# 다중 레이어가 있으면 병합 (v5에서 자주 발생하는 문제)
if(length(Layers(my_wt[["RNA"]])) > 3) {
    cat("다중 레이어 감지됨. 병합 중...\n")
    my_wt[["RNA"]] <- JoinLayers(my_wt[["RNA"]])
}

# 정규화 확인 및 실행
if(!"data" %in% Layers(my_wt[["RNA"]])) {
    cat("정규화 실행 중...\n")
    my_wt <- NormalizeData(my_wt)
}

# ──────────────────────────────────────────────────────────
# 4) 정규화 상태 확인 및 적절한 layer 선택
# ──────────────────────────────────────────────────────────
# 현재 data layer의 값 범위 확인
data_sample <- my_wt[["RNA"]]$data[1:100, 1:5]
data_range <- range(data_sample, na.rm = TRUE)
counts_sample <- my_wt[["RNA"]]$counts[1:100, 1:5]
counts_range <- range(counts_sample, na.rm = TRUE)

cat("Data layer 범위:", round(data_range[1], 3), "~", round(data_range[2], 3), "\n")
cat("Counts layer 범위:", counts_range[1], "~", counts_range[2], "\n")

# 이미 정규화된 데이터인지 확인 (log 스케일이면 보통 0~15 범위)
if (data_range[2] < 20 && data_range[1] >= 0) {
    message("✓ 데이터가 이미 log 정규화되어 있습니다. 기존 data layer 사용.")
} else {
    message("데이터가 정규화되지 않았거나 예상 범위를 벗어남. 확인 필요.")
}

# ──────────────────────────────────────────────────────────
# 5) 유전자 존재 여부 확인
# ──────────────────────────────────────────────────────────
available_genes <- intersect(gene_set, rownames(my_wt))
missing_genes <- setdiff(gene_set, rownames(my_wt))

cat("분석 가능한 유전자 수:", length(available_genes), "/", length(gene_set), "\n")
if (length(missing_genes) > 0) {
    cat("누락된 유전자:", paste(missing_genes, collapse = ", "), "\n")
}

# ──────────────────────────────────────────────────────────
# 6) 평균 발현 계산 (Seurat v5 호환 - 핵심 수정)
# ──────────────────────────────────────────────────────────
if (length(available_genes) == 0) {
    stop("분석 가능한 유전자가 없습니다!")
}

# 클러스터별 평균 발현 계산 (정확한 평균값을 위해 수동 계산)
cat("클러스터별 평균 발현 계산 중...\n")

# 각 클러스터별로 평균 계산
cluster_means <- list()
cluster_names <- unique(my_wt$Myeloid_labels)

for(cluster in cluster_names) {
    # 클러스터별 세포 추출
    cluster_cells <- WhichCells(my_wt, idents = cluster)
    
    if(length(cluster_cells) > 0) {
        # 해당 클러스터 세포들의 평균 발현 계산
        cluster_expr <- rowMeans(GetAssayData(my_wt, slot = "data")[available_genes, cluster_cells, drop = FALSE])
        cluster_means[[cluster]] <- cluster_expr
    }
}

# 데이터프레임으로 변환
avg_mat <- do.call(cbind, cluster_means)
rownames(avg_mat) <- available_genes

# 결과 값 범위 확인 (이제 0-15 범위여야 함)
avg_range <- range(avg_mat, na.rm = TRUE)
cat("평균 발현값 범위:", round(avg_range[1], 3), "~", round(avg_range[2], 3), "\n")

# ──────────────────────────────────────────────────────────
# 7) 데이터 정리 및 시각화
# ──────────────────────────────────────────────────────────
avg_df <- avg_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%
    pivot_longer(cols = -Gene, names_to = "Group", values_to = "Expression") %>%
    mutate(
        Category = case_when(
            Gene %in% wnt_features ~ "Wnt Ligand",
            Gene %in% rct_features ~ "Receptor (rct)",
            Gene %in% tf_features ~ "TF (tf)",
            TRUE ~ "Other"
        ),
        Gene = factor(Gene),
        Expression_rounded = round(Expression, 3)
    )

# 데이터 미리보기
cat("\n상위 발현 유전자:\n")
print(avg_df %>% arrange(desc(Expression)) %>% head(10))

# ──────────────────────────────────────────────────────────
# 8) 클러스터별 barplot 시각화 (개별 그래프들)
# ──────────────────────────────────────────────────────────
# 클러스터별로 개별 barplot 생성
cluster_names <- unique(my_wt$Myeloid_labels)
plot_list <- list()

for(cluster in cluster_names) {
    cluster_data <- avg_df %>% 
        filter(Group == cluster) %>%
        arrange(desc(Expression))
    
    p <- ggplot(cluster_data,
                aes(x = fct_reorder(Gene, Expression),
                    y = Expression,
                    fill = Category)) +
        geom_col(width = 0.8, colour = "black", alpha = 0.8) +
        scale_fill_manual(
            values = c("Wnt Ligand" = "#2ca02c",
                       "Receptor (rct)" = "#4575b4", 
                       "TF (tf)" = "#d73027"),
            name = "Gene Type"
        ) +
        labs(
            title = paste("Wnt Pathway Genes -", cluster),
            x = "Gene",
            y = "Average expression"
        ) +
        theme_few(base_size = 10) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "none",  # 개별 그래프에서는 범례 제거
            plot.title = element_text(size = 11, face = "bold")
        )
    
    plot_list[[cluster]] <- p
}

# patchwork로 모든 그래프 결합
library(patchwork)
combined_plot <- wrap_plots(plot_list, ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")

print(combined_plot)

# ──────────────────────────────────────────────────────────
# 9) 선택사항: 요약 그래프들 (모든 클러스터 통합)
# ──────────────────────────────────────────────────────────
# Wnt Ligand 통합 그래프
wnt_summary <- avg_df %>% filter(Category == "Wnt Ligand")
p_wnt_summary <- ggplot(wnt_summary,
                        aes(x = fct_reorder(Gene, Expression),
                            y = Expression,
                            fill = Category)) +
    geom_col(width = 0.8, colour = "black", alpha = 0.8) +
    facet_wrap(~ Group, scales = "free_y", ncol = 2) +
    scale_fill_manual(
        values = c("Wnt Ligand" = "#0fa0f0"),
        name = "Gene Type"
    ) +
    labs(
        title = "Wnt Ligands Expression by Myeloid Subsets",
        subtitle = paste0("WT cells (n=", ncol(my_wt), " cells)"),
        x = "Gene",
        y = "Average normalized expression"
    ) +
    theme_few(base_size = 10) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.title = element_text(size = 11),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        legend.position = "top",
        strip.text = element_text(size = 9, face = "bold")
    )

# Receptor + TF 통합 그래프
rct_tf_summary <- avg_df %>% filter(Category %in% c("Receptor (rct)", "TF (tf)"))
p_rct_tf_summary <- ggplot(rct_tf_summary,
                           aes(x = fct_reorder(Gene, Expression),
                               y = Expression,
                               fill = Category)) +
    geom_col(width = 0.8, colour = "black", alpha = 0.8) +
    facet_wrap(~ Group, scales = "free_y", ncol = 2) +
    scale_fill_manual(
        values = c("Receptor (rct)" = "#f00fa0", 
                   "TF (tf)" = "#a0f00f"),
        name = "Gene Type"
    ) +
    labs(
        title = "Wnt Receptors & TFs Expression by Myeloid Subsets",
        subtitle = paste0("WT cells (n=", ncol(my_wt), " cells)"),
        x = "Gene",
        y = "Average normalized expression"
    ) +
    theme_few(base_size = 10) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.title = element_text(size = 11),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        legend.position = "top",
        strip.text = element_text(size = 9, face = "bold")
    )

# 통합 그래프들 저장
wnt_summary_filename <- file.path(output_dir, "Wnt_Ligands_All_Clusters_Summary.tiff")
ggsave(filename = wnt_summary_filename,
       plot = p_wnt_summary,
       device = "tiff",
       width = 14,
       height = 10,
       dpi = 600,
       compression = "lzw")

rct_tf_summary_filename <- file.path(output_dir, "Wnt_Receptors_TF_All_Clusters_Summary.tiff")
ggsave(filename = rct_tf_summary_filename,
       plot = p_rct_tf_summary,
       device = "tiff",
       width = 14,
       height = 10,
       dpi = 600,
       compression = "lzw")

cat("통합 요약 그래프 저장 완료:\n")
cat("- Wnt Ligands Summary:", wnt_summary_filename, "\n")
cat("- Receptors + TF Summary:", rct_tf_summary_filename, "\n")

print(p_wnt_summary)
print(p_rct_tf_summary)

# ──────────────────────────────────────────────────────────
# 10) 수치 요약 테이블
# ──────────────────────────────────────────────────────────
# 클러스터별 요약
cluster_summary <- avg_df %>%
    group_by(Group, Category) %>%
    summarise(
        N_genes = n(),
        Mean_expression = round(mean(Expression), 3),
        SD_expression = round(sd(Expression), 3),
        Min_expression = round(min(Expression), 3),
        Max_expression = round(max(Expression), 3),
        .groups = "drop"
    )

cat("\n=== 클러스터별 유전자 그룹 발현 요약 ===\n")
print(cluster_summary)

# 전체 요약
summary_table <- avg_df %>%
    group_by(Category) %>%
    summarise(
        N_genes = n(),
        Mean_expression = round(mean(Expression), 3),
        SD_expression = round(sd(Expression), 3),
        Min_expression = round(min(Expression), 3),
        Max_expression = round(max(Expression), 3),
        .groups = "drop"
    )

cat("\n=== 전체 유전자 그룹별 발현 요약 ===\n")
print(summary_table)

cat("\n✓ 분석 완료!\n")
cat("주요 결과:\n")
cat("- 총", length(available_genes), "개 Wnt pathway 유전자 분석\n")
cat("- 발현 범위:", round(avg_range[1], 3), "~", round(avg_range[2], 3), "\n")
cat("- 가장 높은 발현:", avg_df$Gene[which.max(avg_df$Expression)], 
    "(", round(max(avg_df$Expression), 3), ")\n")





# ──────────────────────────────────────────────────────────
#  S10) 클러스터별 · 카테고리별 개별 그래프 & raw data 저장
# ──────────────────────────────────────────────────────────
# (위에서 avg_df를 이미 계산했다고 가정)
library(ggplot2)
library(dplyr)
library(forcats)

# 1) 출력 디렉터리 설정
output_dir <- "./figure/250610"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2) 클러스터 리스트
cluster_names <- unique(avg_df$Group)

# 3) 색상 맵 (기존과 동일)
col_map <- c(
    "Wnt Ligand"     = "#0fa0f0",
    "Receptor (rct)" = "#f00fa0",
    "TF (tf)"        = "#a0f00f"
)

# 4) 유전자 순서 정의 (사용자 지정 순서)
wnt_order <- c("Wnt1", "Wnt2", "Wnt2b", "Wnt3", "Wnt4", "Wnt5a", "Wnt5b", "Wnt6",
               "Wnt7a", "Wnt7b", "Wnt8a", "Wnt8b", "Wnt9a", "Wnt9b", "Wnt10a", "Wnt10b",
               "Wnt11", "Wnt16")

# 리셉터 순서 (입력 순서대로)
receptor_order <- c('Lrp1','Lrp5','Lrp6',
                    'Fzd1','Fzd2','Fzd3','Fzd4','Fzd5',
                    'Fzd6','Ror1','Ror2','Ryk')

# 전사인자 순서 (입력 순서대로)
tf_order <- c('Tcf3','Tcf4','Tcf7','Lef1')

# 리셉터 + 전사인자 통합 순서
receptor_tf_order <- c(receptor_order, tf_order)

# 5) 클러스터별 반복
for (cl in cluster_names) {
    
    # --- a) Wnt Ligand 데이터 / 플롯 ---
    df_wnt <- avg_df %>% 
        filter(Group == cl, Category == "Wnt Ligand") %>%
        # 사용자 정의 순서로 factor 설정
        mutate(Gene = factor(Gene, levels = wnt_order)) %>%
        arrange(Gene)  # 지정된 순서대로 정렬
    
    # CSV 저장
    write.csv(
        df_wnt,
        file   = file.path(output_dir, paste0(cl, "_WntLigand_expression.csv")),
        row.names = FALSE
    )
    
    # barplot (순서 고정, 글씨색 검정)
    p_w <- ggplot(df_wnt,
                  aes(x = Gene,  # fct_reorder 제거하여 순서 고정
                      y = Expression,
                      fill = Category)) +
        geom_col(colour = "black", width = 0.8) +
        scale_fill_manual(values = col_map["Wnt Ligand"], guide = FALSE) +
        labs(
            title = paste0(cl, " — Wnt Ligands"),
            x     = "Gene",
            y     = "Average expression"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                       colour = "black"),  # 글씨색 검정
            axis.text.y = element_text(colour = "black"),  # Y축 글씨색도 검정
            plot.title  = element_text(face = "bold", size = 14, colour = "black"),
            axis.title  = element_text(colour = "black"),  # 축 제목도 검정
            panel.grid = element_blank()
        )
    
    # TIFF 저장
    ggsave(
        filename = file.path(output_dir, paste0(cl, "_WntLigand_barplot.tiff")),
        plot     = p_w,
        device   = "tiff",
        width    = 4, height = 4, dpi = 600, compression = "lzw"  # width 증가
    )
    
    
    # --- b) Receptor + TF 데이터 / 플롯 ---
    df_rt <- avg_df %>% 
        filter(Group == cl, Category %in% c("Receptor (rct)", "TF (tf)")) %>%
        # 사용자 정의 순서로 factor 설정 (리셉터 먼저, 그 다음 TF)
        mutate(Gene = factor(Gene, levels = receptor_tf_order)) %>%
        arrange(Gene)  # 지정된 순서대로 정렬
    
    # CSV 저장
    write.csv(
        df_rt,
        file   = file.path(output_dir, paste0(cl, "_Receptor_TF_expression.csv")),
        row.names = FALSE
    )
    
    # barplot (순서 고정, 글씨색 검정)
    p_r <- ggplot(df_rt,
                  aes(x = Gene,  # fct_reorder 제거하여 순서 고정
                      y = Expression,
                      fill = Category)) +
        geom_col(colour = "black", width = 0.8) +
        scale_fill_manual(values = col_map[c("Receptor (rct)", "TF (tf)")]) +
        labs(
            title = paste0(cl, " — Receptors & TFs"),
            x     = "Gene",
            y     = "Average expression",
            fill  = "Category"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                       colour = "black"),  # 글씨색 검정
            axis.text.y = element_text(colour = "black"),  # Y축 글씨색도 검정
            plot.title  = element_text(face = "bold", size = 14, colour = "black"),
            axis.title  = element_text(colour = "black"),  # 축 제목도 검정
            legend.text = element_text(colour = "black"),  # 범례 글씨색 검정
            legend.title = element_text(colour = "black"), # 범례 제목 검정
            panel.grid = element_blank()
        )
    
    # TIFF 저장
    ggsave(
        filename = file.path(output_dir, paste0(cl, "_Receptor_TF_barplot.tiff")),
        plot     = p_r,
        device   = "tiff",
        width    = 4, height = 4, dpi = 600, compression = "lzw"  # width 증가
    )
    
}

cat("모든 클러스터별 그래프와 데이터가 저장되었습니다!\n")
cat("저장 위치:", output_dir, "\n")
cat("- Wnt Ligand: 사용자 지정 순서로 정렬 (Wnt1 → Wnt16)\n")
cat("- Receptor: 입력 순서대로 정렬 (Lrp1 → Ryk)\n")
cat("- TF: 입력 순서대로 정렬 (Tcf3 → Lef1)\n")
cat("- 모든 텍스트: 검정색으로 설정\n")
cat("- 그래프 크기: 8x4 인치로 조정\n")


# ──────────────────────────────────────────────────────────
#  S10) 클러스터별 · 카테고리별 개별 그래프 & raw data 저장
# ──────────────────────────────────────────────────────────

# (위에서 avg_df를 이미 계산했다고 가정)

library(ggplot2)
library(dplyr)
library(forcats)

# 1) 출력 디렉터리 설정
output_dir <- "./figure/250610"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2) 클러스터 리스트
cluster_names <- unique(avg_df$Group)

# 3) 색상 맵 (기존과 동일)
col_map <- c(
    "Wnt Ligand"     = "#0fa0f0",
    "Receptor (rct)" = "#f00fa0",
    "TF (tf)"        = "#a0f00f"
)

# 4) 클러스터별 반복
for (cl in cluster_names) {
    
    # --- a) Wnt Ligand 데이터 / 플롯 ---
    df_wnt <- avg_df %>% 
        filter(Group == cl, Category == "Wnt Ligand") %>%
        arrange(desc(Expression))
    
    # CSV 저장
    write.csv(
        df_wnt,
        file   = file.path(output_dir, paste0(cl, "_WntLigand_expression.csv")),
        row.names = FALSE
    )
    
    # barplot
    p_w <- ggplot(df_wnt,
                  aes(x = fct_reorder(Gene, Expression),
                      y = Expression,
                      fill = Category)) +
        geom_col(colour = "black", width = 0.8) +
        scale_fill_manual(values = col_map["Wnt Ligand"], guide = FALSE) +
        labs(
            title = paste0(cl, " — Wnt Ligands"),
            x     = "Gene",
            y     = "Average expression"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(face = "bold", size = 14)
        )

    # TIFF 저장
    ggsave(
        filename = file.path(output_dir, paste0(cl, "_WntLigand_barplot.tiff")),
        plot     = p_w,
        device   = "tiff",
        width    = 6, height = 4, dpi = 600, compression = "lzw"
    )
    
    
    # --- b) Receptor + TF 데이터 / 플롯 ---
    df_rt <- avg_df %>% 
        filter(Group == cl, Category %in% c("Receptor (rct)", "TF (tf)")) %>%
        arrange(Category, desc(Expression))
    
    # CSV 저장
    write.csv(
        df_rt,
        file   = file.path(output_dir, paste0(cl, "_Receptor_TF_expression.csv")),
        row.names = FALSE
    )
    
    # barplot
    p_r <- ggplot(df_rt,
                  aes(x = fct_reorder(Gene, Expression),
                      y = Expression,
                      fill = Category)) +
        geom_col(colour = "black", width = 0.8, position = position_dodge(width = 0.9)) +
        scale_fill_manual(values = col_map[c("Receptor (rct)", "TF (tf)")]) +
        labs(
            title = paste0(cl, " — Receptors & TFs"),
            x     = "Gene",
            y     = "Average expression",
            fill  = "Category"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(face = "bold", size = 14)
        )
    
    # TIFF 저장
    ggsave(
        filename = file.path(output_dir, paste0(cl, "_Receptor_TF_barplot.tiff")),
        plot     = p_r,
        device   = "tiff",
        width    = 6, height = 4, dpi = 600, compression = "lzw"
    )
    
}


## ───────────────────────────────────────────────────────────────────────
## DotPlot : 원하는 y-축(세포 타입) 순서 지정 버전
##   · Cancer cell 계열 → CAF → Endothelial → Myeloid/Innate → T/NK → B
##   · 기존 colour palette · feature 목록 그대로 사용
## ───────────────────────────────────────────────────────────────────────

library(Seurat);  library(ggplot2);  library(RColorBrewer)

## 1) y-축(= Idents) 순서 정의  -----------------------------------------
desired_order <- c(
    ## Cancer cells (맨 위로)
    "Gli3hi Cancer cell", "Tffhi Cancer cell",
    "Mki67intTop2aint Cancer cell", "Mki67hiTop2ahi Cancer cell",
    "Idohi cancer cell", "Epi/Cancer cell",
    
    ## T / NK
    "NK cell", "CD4 T cell", "Mki67hiTop2Ahi CTL", "CTL",
    
    ## CAF
    "Pdgfahi CAF", "Dcnhi CAF", "Cdh11hi CAF", "Apoehi CAF",
    
    ## Myeloid / Innate
    "Mki67hiTop2AhiBatf3hi DC", "Batf3hi DC", "Mono/Mac", "Neutrophil",
    
    ## B-cell
    "Mki67hiTop2Ahi B cell", "Epcamhi B cell", "B cell",
    
    ## Endothelial (맨 아래)
    "Endothelial cell"
)

## 2) Idents 순서 적용  --------------------------------------------------
porcn.combined.harmony$NH_labels <-
    factor(porcn.combined.harmony$NH_labels, levels = desired_order)

Idents(porcn.combined.harmony) <- "NH_labels"   # DotPlot 은 Idents 순서를 그대로 그림

## 3) DotPlot  -----------------------------------------------------------
total_features <-  unique(c(
    'Ptprc', 'Krt19', 'Cdh1', 'Msln', 
    'Mki67','Tff1','Tff2','Tff3','Gli3','Top2a', 'Ido1', 'Ido2','Epi', # cancer cell
    'Nkg7', 'Klrg1',
    'Cd3d','Cd3e', 'Cd4', 'Cd8a', 'Sell',
    'Ccr7', 'Il7r', 'Cd28', 'Il2ra', 'Il1rl1',
    'Gata3', 'Trdc', 'Col1a2', 'Pdpn', 
    'Pdgfa', 'Dcn','Cdh11', 'Apoe', #CAF
    'Cd68', 'Adgre1',
    'Itgam', 'Cd14', 'Mrc1', 'H2-Eb1', 'Batf3','S100a8', #DC
    'Batf3',
    'Cd79a', 'Cd19', 'Ms4a1', 'Ly6c2', 'Ly6g',
    'Mrc1', 'H2-Eb1', 'Krt19', 'Cdh1', 
    'Epcam',
    'Try4', 'Pecam1','Cdh5'
))

dot <- DotPlot(
    object    = porcn.combined.harmony,
    features  = total_features,
    cols      = "RdYlBu",
    dot.scale = 10
) +
    RotatedAxis() +
    theme_bw(base_size = 25) +
    theme(
        axis.text.x    = element_text(angle = 90, hjust = .5, vjust = .5,
                                      colour = "black", size = 20),
        axis.text.y    = element_text(colour = "black", size = 20),
        axis.ticks     = element_blank(),
        legend.title   = element_text(size = 18, face = "bold"),
        legend.text    = element_text(size = 18),
        panel.border   = element_rect(colour = "black", fill = NA, size = 1.3)
    ) +
    labs(x = NULL, y = NULL)

ggsave("./figure/250610_annotation_gene_dotplot_reordered.tiff",
       plot = dot,
       width = 20, height = 10, dpi = 600, compression = "lzw")
