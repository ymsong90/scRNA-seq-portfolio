# scRNA-seq Analysis Pipeline - Quick Start Guide

## 📋 Overview

이 분석 파이프라인은 PORCN KO vs WT 비교 연구를 위한 완전한 single-cell RNA-seq 분석 워크플로우입니다.

## 🚀 Quick Start

### 1. 환경 설정

```r
# 필수 패키지 설치
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "tidyr", "tibble"))

# Bioconductor 패키지
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"))

# Harmony (batch correction)
devtools::install_github("immunogenomics/harmony")
```

### 2. 데이터 준비

```bash
# 프로젝트 디렉토리 구조
scRNA-seq-portfolio/
├── data/
│   ├── wt/filtered_feature_bc_matrix/
│   └── ko/filtered_feature_bc_matrix/
├── results/  # 자동 생성됨
└── [분석 스크립트들...]
```

### 3. 순차적 실행

```r
# Step 01: 데이터 로딩 및 QC
source("01_preprocessing/01_data_loading_QC.R")

# Step 02: Normalization 및 Harmony integration
source("02_integration/02_harmony_integration.R")

# Step 03: Clustering 및 annotation
source("03_clustering/03_clustering_annotation.R")

# Step 04: Cell proportion 분석
source("04_visualization/04_cell_proportion_analysis.R")

# Step 05: Differential expression
source("05_differential_expression/05a_myeloid_DEG_volcano.R")
source("05_differential_expression/05b_CD8T_DEG_analysis.R")

# Step 06: Functional enrichment
source("06_functional_enrichment/06a_myeloid_GO_analysis.R")
```

## 📁 파일 구조 상세

### 분석 파이프라인
```
01_preprocessing/
├── 01_data_loading_QC.R          # 데이터 로딩 및 QC
└── README.md                      # 상세 설명

02_integration/
├── 02_harmony_integration.R      # Batch correction
└── README.md

03_clustering/
├── 03_clustering_annotation.R    # Cell type annotation
└── README.md

04_visualization/
├── 04_cell_proportion_analysis.R # 세포 비율 분석
└── README.md

05_differential_expression/
├── 05a_myeloid_DEG_volcano.R     # Myeloid DEG 분석
├── 05b_CD8T_DEG_analysis.R       # CD8 T cell DEG 분석
└── README.md

06_functional_enrichment/
├── 06a_myeloid_GO_analysis.R     # Myeloid GO 분석
├── 06b_monocyte_macrophage_GO.R  # Mono/Mac GO 분석
├── 06c_wnt_pathway_analysis.R    # Wnt pathway 분석
└── README.md

utils/
├── plotting_functions.R           # 공통 plotting 함수
├── deg_functions.R                # DEG 분석 함수
└── go_functions.R                 # GO 분석 함수
```

### 결과 파일
```
results/
├── 01_QC/
│   ├── QC_violin_plots_before_filtering.png
│   └── QC_violin_plots_after_filtering.png
├── 02_integration/
│   ├── UMAP_harmony_clusters.png
│   └── harmony_convergence.png
├── 03_clustering/
│   ├── UMAP_annotated.png
│   ├── DotPlot_markers.tiff
│   └── cluster_markers_all.csv
├── 04_visualization/
│   └── cell_proportion_plots/
├── 05_differential_expression/
│   ├── DEG/
│   ├── volcano/
│   └── ...
└── 06_functional_enrichment/
    ├── GO/
    └── plots/
```

## 🔧 커스터마이징

### QC 파라미터 조정

```r
# 01_preprocessing/01_data_loading_QC.R 수정
QC_PARAMS <- list(
    min_features = 200,      # 최소 유전자 수
    max_features = 8000,     # 최대 유전자 수 (doublet 제거)
    max_mt_pct   = 20        # 최대 미토콘드리아 비율
)
```

### DEG 기준 조정

```r
# 05_differential_expression/*.R 수정
PARAMS <- list(
    deg = list(
        logfc = 0.25,        # Log fold-change 기준값
        minpct = 0.10,       # 최소 detection rate
        padj = 0.05          # Adjusted p-value cutoff
    )
)
```

### GO 분석 파라미터

```r
# 06_functional_enrichment/*.R 수정
go_params <- list(
    pval = 0.10,            # P-value cutoff
    qval = 0.20,            # Q-value cutoff
    ont = "BP",             # Ontology: BP, MF, CC
    showN = 10              # Top N pathways to show
)
```

## 📊 예상 분석 시간

| Step | 예상 시간 | 메모리 사용량 |
|------|----------|-------------|
| 01: QC | 5-10분 | ~4GB |
| 02: Integration | 10-20분 | ~8GB |
| 03: Clustering | 15-30분 | ~10GB |
| 04: Visualization | 5-10분 | ~4GB |
| 05: DEG Analysis | 20-40분 | ~8GB |
| 06: GO Analysis | 10-20분 | ~6GB |

*예상 시간은 ~30,000 cells 기준

## 💡 Tips & Best Practices

### 1. 메모리 관리
```r
# 분석 중간중간 메모리 정리
rm(unnecessary_objects)
gc()

# 큰 객체 저장 후 제거
save(large_object, file = "checkpoint.RData")
rm(large_object)
```

### 2. 체크포인트 활용
각 단계마다 RData 파일이 저장되므로, 중간부터 재시작 가능:
```r
# Step 03부터 시작하려면
load("./data/porcn.combined.harmony.RData")
source("03_clustering/03_clustering_annotation.R")
```

### 3. 병렬 처리
```r
# DEG 분석 속도 향상
library(future)
plan(multisession, workers = 8)  # CPU 코어 수에 맞게 조정
```

## 🐛 Troubleshooting

### 문제: "Error: Cannot find X in object"
**해결**: 이전 단계의 결과 파일이 제대로 저장되었는지 확인
```r
# 저장된 객체 확인
load("./data/porcn.combined.harmony.RData")
names(porcn.combined.harmony@meta.data)
```

### 문제: 메모리 부족 에러
**해결**: 
1. 불필요한 R 프로세스 종료
2. 분석을 더 작은 서브셋으로 분할
3. 서버에서 실행 (권장: 32GB+ RAM)

### 문제: Harmony가 수렴하지 않음
**해결**:
```r
# 더 많은 iteration 허용
porcn.combined.harmony <- RunHarmony(
    porcn.combined,
    "ID",
    max.iter.harmony = 50  # 기본값 10에서 증가
)
```

## 📧 면접 대비 포인트

### 코드 설명 준비
1. **QC 필터링 기준 선택 이유**
2. **Harmony vs CCA/RPCA 선택 이유**
3. **Wilcoxon test 사용 이유**
4. **Multiple testing correction 방법**

### 결과 해석 준비
1. **주요 cell type별 특징**
2. **KO vs WT 차이의 생물학적 의미**
3. **Unexpected findings 및 해석**
4. **Technical limitations & future directions**

### 실행 데모 준비
```r
# 빠른 데모용 (5분 이내)
# 1. 데이터 로딩
load("./data/porcn.combined.harmony_annotated.RData")

# 2. UMAP 시각화
DimPlot(porcn.combined.harmony, group.by = "NH_labels", label = TRUE)

# 3. 특정 마커 확인
FeaturePlot(porcn.combined.harmony, features = c("Cd68", "Cd3e"))

# 4. DEG 결과 확인
deg_results <- read.csv("./results/05_differential_expression/DEG/Myeloid_DEG_ko_vs_wt.csv")
head(deg_results)
```

## 📚 추가 리소스

- **Seurat Tutorials**: https://satijalab.org/seurat/
- **Harmony Documentation**: https://github.com/immunogenomics/harmony
- **clusterProfiler Book**: https://yulab-smu.top/biomedical-knowledge-mining-book/

## 🤝 GitHub Repository 구성

```bash
# .gitignore 설정 (권장)
data/                    # 원본 데이터는 제외
*.RData                  # 큰 R 객체 제외
*.rds
results/                 # 결과 파일 제외 (또는 대표 예시만 포함)

# README.md에 포함할 내용
- 프로젝트 개요
- 분석 파이프라인 설명
- 주요 결과 요약
- 사용 방법
- 예시 output 몇 개
```

---

**Note**: 이 가이드는 면접 준비용으로 작성되었습니다. 실제 분석 시 데이터 특성에 맞게 파라미터를 조정하세요.
