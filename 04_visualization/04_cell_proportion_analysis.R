# =============================================================================
# 필요한 패키지 로드
# =============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)    # percent_format() 사용
library(forcats)   # fct_relevel() 등 factor 순서 재정렬
library(dplyr)

# =============================================================================
# 0) 사전 준비
#   - porcn.combined.harmony 객체가 메모리에 있어야 합니다.
#   - 메타데이터 컬럼 중, “NH_labels” (예: “Classical Monocyte” 등)와
#     “ID” (예: “porcn_wt” / “porcn_ko”) 가 반드시 존재해야 합니다.
# =============================================================================

# =============================================================================
# 1) 메타데이터 추출 및 전처리
# =============================================================================
meta_all <- porcn.combined.harmony@meta.data %>%
    as_tibble(rownames = "cell") %>%
    # 필요한 칼럼만 골라서 이름 바꾸기
    dplyr::select(
        cell, 
        cluster   = NH_labels,   # 세포 타입 라벨
        condition = ID           # WT vs KO 구분
    ) %>%
    # factor 레벨 지정
    mutate(
        # condition 컬럼: “porcn_wt” → “WT”, “porcn_ko” → “KO”
        condition = case_when(
            condition == "porcn_wt" ~ "WT",
            condition == "porcn_ko" ~ "KO",
            TRUE                    ~ as.character(condition)
        ) %>% factor(levels = c("WT", "KO")),
        
        # cluster(NH_labels) 순서는 원하는 순서가 있으면 fct_relevel() 로 지정
        # 여기서는 예시로 알파벳 순서 그대로 사용하거나, 필요 시 순서를 고쳐주세요.
        cluster = factor(cluster) 
    )

# 확인: 어떤 라벨들이 있는지, WT/KO가 어떻게 들어 있는지
table(meta_all$cluster, meta_all$condition)

# =============================================================================
# 2) “세포 타입별 WT vs KO 비율” 바 플롯
#    → x축: cluster (NH_labels), fill: condition (WT/KO), 위치는 “fill” 사용하여 비율 표시
# =============================================================================
library(scales)

p1 <- ggplot(meta_all, aes(x = cluster, fill = condition)) +
    geom_bar(position = "fill", color = "grey30", width = 0.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(
        values = c("WT" = "#e81666", "KO" = "#16e898"),
        labels = c("WT", "KO")
    ) +
    theme_minimal(base_size = 20) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20, colour = 'black'),
        axis.ticks.x       = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
    ) 

p1_1 <- ggplot(meta_all, aes(x = cluster, fill = condition)) +
    geom_bar(position = "fill", color = "grey30", width = 0.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(
        values = c("WT" = "#e81666", "KO" = "#16e898"),
        labels = c("WT", "KO")
    ) +
    theme_minimal(base_size = 20) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
    ) + NoLegend()

# 플롯 출력
print(p1_1)

# -----------------------------------------------------------------------------
# 2-1) 위 플롯의 비율 정보를 CSV로도 저장하고 싶다면
#     → p1$data 에서 count와 비율 계산
# -----------------------------------------------------------------------------
df1 <- as.data.frame(p1$data)

cluster_cond_props1 <- df1 %>%
    group_by(cluster, condition) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()

# CSV 저장
write.csv(
    cluster_cond_props1,
    file      = "./data/250604_porcn_all_NHlabels_WT_vs_KO_props_by_cluster.csv",
    row.names = FALSE,
    quote     = TRUE
)

# -----------------------------------------------------------------------------
# 3) “WT vs KO 내에서 전체 세포 타입(NH_labels) 분포” 스택형 바 플롯
# -----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)

# --- 1) meta_all 객체 준비 ----------------------------------------------------
# porcn.combined.harmony@meta.data 에는 다음 두 칼럼이 있어야 합니다:
#   • NH_labels : 세포 타입 (factor)
#   • ID        : “porcn_wt” 또는 “porcn_ko” (character 또는 factor)

meta_all <- porcn.combined.harmony@meta.data %>%
    as_tibble(rownames = "cell") %>%
    dplyr::select(
        cell,
        cluster   = NH_labels,   # 전체 세포 타입 라벨
        condition = ID           # “porcn_wt” vs “porcn_ko”
    ) %>%
    mutate(
        # “porcn_wt” → “WT”, “porcn_ko” → “KO”로 바꾸고 factor 순서 지정
        condition = case_when(
            condition == "porcn_wt" ~ "WT",
            condition == "porcn_ko" ~ "KO",
            TRUE                    ~ as.character(condition)
        ) %>% factor(levels = c("WT", "KO")),
        
        # cluster는 NH_labels에 이미 모든 레벨이 factor로 존재하므로,
        # 필요하다면 levels 순서를 fct_relevel()로 재지정할 수 있습니다.
        # 여기서는 원래 순서 그대로 유지합니다.
        cluster = factor(cluster)
    )

# --- 2) 색상 벡터 정의 --------------------------------------------------------
# porcn.combined.harmony에 있는 NH_labels 레벨에 맞춰 colors 벡터를 작성합니다.
# (아래는 예시로 제공된 celltype.cols를 그대로 사용)
celltype.cols <- c(
    # ----- B-cell 계열 -----
    "B cell"                    = "#808000",
    "Epcamhi B cell"            = "#B1D337",
    "Mki67hiTop2Ahi B cell"     = "#21A04A",
    
    # ----- Neutrophil / DC / Macrophage -----
    "Neutrophil"                  = "#778899",
    "Batf3hi DC"                = "#F17D97",
    "Mki67hiTop2AhiBatf3hi DC"  = "#A7213A",
    "Mono/Mac"                  = "#ED4169",
    
    # ----- Cancer-cell 계열 -----
    "Epi/Cancer cell"           = "#0B6A42",
    "Idohi cancer cell"         = "#CDE9C8",
    "Mki67intTop2aint Cancer cell" = "#0B6B6A",
    "Tffhi Cancer cell"         = "#22A1B1",
    "Mki67hiTop2ahi Cancer cell"= "#73CCD4",
    "Gli3hi Cancer cell"        = "#C2E5EB",
    
    # ----- CAF 계열 -----
    "Apoehi CAF"                = "#F65C4B",
    "Cdh11hi CAF"               = "#F68F80",
    "Dcnhi CAF"                 = "#EA6D22",
    "Pdgfahi CAF"               = "#F58927",
    
    # ----- Endothelial -----
    "Endothelial cell"          = "#F0B817",
    
    # ----- CTL / CD4 T -----
    "CTL"                       = "#C2E9F9",
    "Mki67hiTop2Ahi CTL"        = "#82D5F7",
    "CD4 T cell"                = "#0E7DC2",
    
    # ----- NK -----
    "NK cell"                   = "#8F7BB7"
)

# (※ 만약 porcn.combined.harmony$NH_labels에 위 목록 외 다른 레벨이 있다면,
#     반드시 모든 레벨에 대응하는 색상을 추가/수정해 주세요.)

# --- 3) 스택형 바 플롯 그리기 --------------------------------------------------
# 1) cluster 순서 설정 -----------------------------------------
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

# cluster 순서를 원하는 대로 설정
meta_all$cluster <- factor(meta_all$cluster, levels = desired_order)

# 2) ggplot 코드 -----------------------------------------------
p_all <- ggplot(meta_all, aes(x = condition, fill = cluster)) +
    
    # WT/KO 각각을 0–1 비율로 스케일링
    geom_bar(position = "fill", width = 0.6, color = "grey30") +
    
    # y축을 백분율로 표시
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    
    # cluster별 색상 매핑 (celltype.cols 사용)
    scale_fill_manual(values = celltype.cols) +
    
    # 축/범례 레이블
    labs(
        x    = "Condition",
        y    = "Proportion (%)",
        fill = "Cell Type (NH_labels)"
    ) +
    
    # 테마 설정: 불필요한 그리드 제거, 축 텍스트만 남김
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x       = element_text(face = "bold", size = 14),
        axis.title.y       = element_text(face = "bold", size = 14),
        axis.text.x        = element_text(size = 12, colour = "black"),
        axis.text.y        = element_text(size = 12, colour = "black"),
        axis.ticks.x       = element_blank()
    )

# 3) 그래프 출력
print(p_all)

p_all <- ggplot(meta_all, aes(x = condition, fill = cluster)) +
    
    # WT/KO 각각을 0–1 비율로 스케일링
    geom_bar(position = "fill", width = 0.6, color = "grey30") +
    
    # y축을 백분율로 표시
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    
    # cluster별 색상 매핑 (celltype.cols 사용)
    scale_fill_manual(values = celltype.cols) +
    
    # 축/범례 레이블
    labs(
        x    = "Condition",
        y    = "Proportion (%)",
        fill = "Cell Type (NH_labels)"
    ) +
    
    # 테마 설정: 불필요한 그리드 제거, 축 텍스트만 남김
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x       = element_blank(),
        axis.title.y       = element_blank(),
        axis.text.x        = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.x       = element_blank()
    ) + NoLegend()

# 플롯 확인
print(p_all)

# --- 4) “WT vs KO” 내에서 각 Cell Type 비율 정보 CSV로 저장 ----------------------
df_all <- as.data.frame(p_all$data)

cluster_cond_props_all <- df_all %>%
    group_by(condition, cluster) %>%
    summarise(
        count = n(),
        .groups = "drop"
    ) %>%
    group_by(condition) %>%
    mutate(
        prop = count / sum(count)
    ) %>%
    ungroup()

# CSV 파일로 저장
write.csv(
    cluster_cond_props_all,
    file      = "./figure/250604_porcn_all_NHlabels_props_WT_KO.csv",
    row.names = FALSE,
    quote     = TRUE
)



# =============================================================================
# 4) 최종 출력물 저장 (TIFF, 600 DPI)
# =============================================================================

# 4-1) 세포 타입별 WT vs KO 비율 바 플롯 저장
ggsave(
    filename = "./figure/250604_porcn_all_NHlabels_WT_vs_KO_by_cluster_forpaper.tiff",
    plot     = p1_1,
    device   = "tiff",
    dpi      = 600,
    width    = 15,   # 필요에 따라 조절 (inch 단위)
    height   = 10,
    units    = "in"
)

# 4-2) WT vs KO 내 세포 타입 분포 스택형 바 플롯 저장
ggsave(
    filename = "./figure/250602_porcn_all_NHlabels_stackbar_by_condition.tiff",
    plot     = p2,
    device   = "tiff",
    dpi      = 600,
    width    = 4.5, # 필요에 따라 조절
    height   = 5,
    units    = "in"
)

# 4-3) 
ggsave(
    filename = "./figure/250604_porcn_all_NHlabels_stackbar_by_condition_nolabel.tiff",
    plot     = p_myeloid,
    device   = "tiff",
    dpi      = 600,
    width    = 3, # 필요에 따라 조절
    height   = 10,
    units    = "in"
)

# 4-3) 
ggsave(
    filename = "./figure/250613/250604_porcn_all_NHlabels_stackbar_by_condition_nolabel.tiff",
    plot     = p_all,
    device   = "tiff",
    dpi      = 600,
    width    = 8, # 필요에 따라 조절
    height   = 10,
    units    = "in"
)

ggsave(
    filename = "./figure/250604_porcn_all_NHlabels_stackbar_by_condition.tiff",
    plot     = p_all,
    device   = "tiff",
    dpi      = 600,
    width    = 8, # 필요에 따라 조절
    height   = 10,
    units    = "in"
)


### 250613_bargraph 수정, rawdata 추출

# 1) 각 그룹별 비율 계산 --------------------------------------
# `meta_all` 데이터에서 각 조건(condition)별로 cluster에 대한 비율을 계산
df_summary <- meta_all %>%
    group_by(condition, cluster) %>%
    tally() %>%
    group_by(condition) %>%
    mutate(proportion = n / sum(n))  # 비율 계산

# 2) CSV 파일로 저장 -----------------------------------------
# 비율을 포함한 데이터 프레임을 CSV 파일로 저장
write.csv(df_summary, "./figure/250613/cell_type_proportions.csv", row.names = FALSE)

# CSV 파일을 확인할 수 있게 출력 (옵션)
head(df_summary)
