# Single-Cell RNA-seq Analysis Portfolio

마우스 PORCN KO vs WT 비교 연구 및 Human 췌장암 세포 간 상호작용 분석을 포함한 포괄적인 scRNA-seq 분석 파이프라인입니다.

## 📊 프로젝트 개요

이 레포지토리는 두 가지 독립적인 scRNA-seq 프로젝트의 분석 코드를 포함합니다:

### 메인 프로젝트 (Steps 01-06): Mouse PORCN KO 분석
PORCN knockout 마우스 모델에서 면역 미세환경의 변화를 규명하는 comprehensive scRNA-seq 분석입니다.

**주요 발견:**
- 20개 이상의 distinct cell population 프로파일링
- Myeloid cell subset 특성화 및 functional annotation
- 조건별 differential gene expression 분석
- PORCN 결핍에 의한 생물학적 경로 변화 규명

### 추가 프로젝트 (Step 07): Human Pancreatic Cancer CellChat
Human 췌장암 scRNA-seq 데이터에서 WNT 신호전달을 통한 세포 간 상호작용을 분석합니다.

**주요 내용:**
- Cell-cell communication network 분석
- WNT pathway-mediated interaction 규명
- Myeloid-Epithelial cell crosstalk 특성화

## 🔬 기술 스택

**언어:** R (version ≥ 4.0)

**핵심 패키지:**
- `Seurat` (v5+): Single-cell 분석 프레임워크
- `Harmony`: Batch effect correction
- `CellChat`: Cell-cell communication 분석
- `clusterProfiler`: GO/KEGG enrichment 분석
- `ggplot2`, `patchwork`: Visualization

**분석 규모:**
- ~20,000-50,000 cells
- 2 conditions (WT vs KO)
- 20+ cell type annotations

## 📁 레포지토리 구조

```
scRNA-seq-portfolio/
├── README.md                           # 이 파일
│
├── 01_preprocessing/                   # Mouse: 데이터 로딩 및 QC
│   ├── 01_data_loading_QC.R
│   └── README.md
│
├── 02_integration/                     # Mouse: Batch correction
│   ├── 02_harmony_integration.R
│   └── README.md
│
├── 03_clustering/                      # Mouse: Clustering 및 annotation
│   ├── 03_clustering_annotation.R
│   └── README.md
│
├── 04_visualization/                   # Mouse: 세포 비율 분석
│   ├── 04_cell_proportion_analysis.R
│   └── README.md
│
├── 05_differential_expression/         # Mouse: DEG 분석
│   ├── 05a_myeloid_DEG_volcano.R
│   ├── 05b_CD8T_DEG_analysis.R
│   └── README.md
│
├── 06_functional_enrichment/           # Mouse: GO/KEGG 분석
│   ├── 06a_myeloid_GO_analysis.R
│   ├── 06b_monocyte_macrophage_GO.R
│   ├── 06c_wnt_pathway_analysis.R
│   └── README.md
│
├── 07_cell_communication/              # Human: CellChat 분석
│   ├── 07_cellchat_analysis_human.R
│   └── README.md
│
└── utils/                              # 공통 함수들
    └── plotting_functions.R
```

## 🚀 분석 파이프라인

### Step 01: 데이터 전처리 및 QC
**Dataset:** Mouse PORCN KO vs WT

```r
# 10X Genomics 데이터 로딩
# Quality control (nFeature, nCount, percent.mt)
# 필터링: 200 < nFeature < 8000, percent.mt < 20%
```

**핵심 포인트:**
- Mouse mitochondrial gene pattern 이해 ("mt-" lowercase)
- QC threshold 설정의 생물학적 근거

**Output:** Filtered Seurat object

---

### Step 02: Normalization 및 Integration
**Dataset:** Mouse PORCN KO vs WT

```r
# LogNormalize (scale factor: 10,000)
# Highly variable features 선정 (top 2,000)
# Harmony를 이용한 batch correction (WT/KO)
```

**핵심 포인트:**
- PCA를 위한 scaling의 중요성
- Harmony vs 다른 integration 방법 비교

**Output:** Batch-corrected integrated object

---

### Step 03: Clustering 및 Cell Type Annotation
**Dataset:** Mouse PORCN KO vs WT

```r
# PCA (50 PCs 계산, 1:30 사용)
# UMAP/t-SNE visualization
# Graph-based clustering (Louvain algorithm)
# Canonical marker를 이용한 cell type annotation
```

**주요 Cell Type:**
- **Cancer cells** (6 subtypes)
- **T cells**: CD4 T, CD8 T (CTL subtypes), NK
- **Myeloid cells**: Monocytes, Macrophages, DC, Neutrophils
- **B cells** (3 subtypes)
- **Stromal cells**: CAF (4 subtypes), Endothelial

**특수 분석:**
- Myeloid subset의 high-resolution reclustering
- 5개 distinct monocyte/macrophage population 규명

**Output:** Annotated Seurat object (NH_labels, Complete_Labels)

---

### Step 04: Cell Type Proportion 분석
**Dataset:** Mouse PORCN KO vs WT

```r
# 세포 유형별 분포 계산
# WT vs KO 비율 비교
# 통계적 유의성 검정 (Chi-square)
```

**Visualization:**
- Stacked bar plot
- Log2 fold change plot

**Output:** 비율 표, 통계 결과, 시각화

---

### Step 05: Differential Expression Analysis
**Dataset:** Mouse PORCN KO vs WT

```r
# Cluster-specific KO vs WT 비교
# FindMarkers (Wilcoxon rank-sum test)
# Filtering: |log2FC| > 0.25, p.adj < 0.05, min.pct > 0.10
```

**분석 대상:**
- **05a:** Myeloid cell 5개 subcluster별 DEG
- **05b:** CD8 T cell DEG (KO vs WT, CD8 vs Rest)

**핵심 포인트:**
- Non-parametric test 선택 이유
- Multiple testing correction (Bonferroni)
- Volcano plot을 통한 효과적인 visualization

**Output:** DEG tables, volcano plots

---

### Step 06: Functional Enrichment Analysis
**Dataset:** Mouse PORCN KO vs WT

```r
# GO Biological Process enrichment
# Gene symbol → Entrez ID 변환
# clusterProfiler를 이용한 enrichment
```

**분석 내용:**
- **06a:** Myeloid cell 전체 GO enrichment
- **06b:** Monocyte/Macrophage 특이적 GO
- **06c:** WNT pathway gene 발현 분석

**핵심 포인트:**
- Entrez ID 변환의 필요성
- Multiple ontology 활용 (BP, MF, CC)

**Output:** GO tables, enrichment plots, WNT gene heatmap

---

### Step 07: Cell-Cell Communication Analysis (CellChat)
**⚠️ 별도 프로젝트:** Human Pancreatic Cancer scRNA-seq

```r
# Human ligand-receptor database 사용
# WNT signaling pathway 중심 분석
# Myeloid → Epithelial cell interaction
```

**핵심 포인트:**
- **Species-specific 처리:**
  - Human: Gene name uppercase + Ensembl version 제거
  - Mouse: Gene name 그대로 사용
- **Database 선택:** CellChatDB.human vs CellChatDB.mouse
- **Parameter 조정:** min.cells threshold의 의미

**분석 내용:**
- WNT pathway-mediated cell-cell interaction
- Myeloid cell의 cancer cell behavior 조절
- Ligand-receptor pair 규명

**Output:** CellChat object, network visualizations, interaction tables

**생물학적 의의:**
- Tumor-promoting microenvironment 형성 메커니즘
- 치료 타겟 후보 규명

---

## 💡 코드 품질 특징

### 가독성 (Readability)
- **모듈화된 구조:** 각 분석 단계를 독립적인 스크립트로 분리
- **명확한 naming convention:** 함수명, 변수명에서 의도 명확히 전달
- **간결한 주석:** 중요한 decision point에만 `NOTE:` 주석 추가

### 재현성 (Reproducibility)
- **Fixed parameters:** 모든 threshold와 cutoff 명시
- **명시적 random seed:** 필요시 설정 (현재 코드에는 생략)
- **Input/output 명확화:** 각 스크립트의 의존성 명시

### 유지보수성 (Maintainability)
- **Parameter list 분리:** 수정이 쉬운 구조
- **에러 처리:** 존재하지 않는 gene/cluster 체크
- **Output 정리:** 체계적인 디렉토리 구조

### 예제 - 깔끔한 코드 스타일
```r
################################################################################
# Step 03: Clustering and Cell Type Annotation
#
# Purpose: Identify marker genes and annotate cell types
# Dataset: Mouse PORCN KO vs WT
################################################################################

library(Seurat)
library(dplyr)

# NOTE: Join layers before FindAllMarkers (required for Seurat v5)
seu <- JoinLayers(object = seu)

# Find differentially expressed genes
markers <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)
```

**특징:**
- ✅ 불필요한 `cat()` 출력 제거
- ✅ `suppressPackageStartupMessages()` 제거
- ✅ 중요한 부분에만 `NOTE:` 주석
- ✅ 명확한 헤더와 구조

---

## 📊 주요 분석 결과

### Myeloid Cell Characterization
**5개 distinct population 규명:**
1. Classical Monocyte
2. Activated Monocyte
3. Trem2+ Macrophage
4. Selenop+ Macrophage
5. Hexb+ Macrophage

### PORCN 결핍 효과
**면역 세포 구성 변화:**
- Myeloid cell 비율 변화
- T cell activation state shift

**전사체 reprogramming:**
- WNT pathway gene 유의한 변화
- Immune response gene dysregulation

### Cell-Cell Communication (Human Pancreatic Cancer)
**WNT signaling:**
- Myeloid cell → Epithelial cell WNT ligand 분비
- FZD receptor family 발현 패턴
- Tumor-promoting microenvironment 형성

---

## 💻 사용법

### Prerequisites
```r
# Core packages
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"))

# Integration
install.packages("harmony")

# Cell-cell communication
devtools::install_github("sqjin/CellChat")
```

### 전체 파이프라인 실행
```r
# Mouse PORCN KO Analysis (Steps 01-06)
source("01_preprocessing/01_data_loading_QC.R")
source("02_integration/02_harmony_integration.R")
source("03_clustering/03_clustering_annotation.R")
source("04_visualization/04_cell_proportion_analysis.R")
source("05_differential_expression/05a_myeloid_DEG_volcano.R")
source("05_differential_expression/05b_CD8T_DEG_analysis.R")
source("06_functional_enrichment/06a_myeloid_GO_analysis.R")

# Human Pancreatic Cancer CellChat (Step 07)
# NOTE: 별도 데이터 필요
source("07_cell_communication/07_cellchat_analysis_human.R")
```

### Input Data 구조 (Mouse)
```
data/
├── wt/
│   └── filtered_feature_bc_matrix/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
└── ko/
    └── filtered_feature_bc_matrix/
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
```

---

## 📈 예제 Output

### UMAP Visualization
- Cell type annotation이 표시된 UMAP projection
- Condition별 split view (WT vs KO)

### Cell Proportion Analysis
- Stacked bar plot: 조건별 세포 유형 비율
- Log2 fold change plot: 변화량 시각화

### Volcano Plots
- Cluster별 KO vs WT DEG 시각화
- Top significant genes 자동 labeling

### GO Enrichment Plots
- Upregulated/downregulated genes의 enriched pathway
- Dot plot, bar plot 형태

### CellChat Network (Human)
- Circle plot: 전체 WNT 네트워크
- Chord diagram: Directional communication

---

## 🔧 주요 분석 Parameters

### Quality Control (Step 01)
```r
min_features = 200      # Cell당 최소 gene 수
max_features = 8000     # Cell당 최대 gene 수
max_mt_pct   = 20       # 최대 mitochondrial %
```

### Integration (Step 02)
```r
n_variable_features = 2000  # Highly variable genes
n_pcs = 50                  # PCA components
use_pcs = 1:30              # Clustering에 사용할 PCs
```

### Differential Expression (Step 05)
```r
logfc_threshold = 0.25  # ~1.2배 fold change
min_pct = 0.10          # 최소 10% 세포에서 발현
padj_cutoff = 0.05      # Adjusted p-value
test = "wilcox"         # Wilcoxon rank-sum test
```

### GO Enrichment (Step 06)
```r
ont = "BP"              # Biological Process
pvalueCutoff = 0.05
qvalueCutoff = 0.05
OrgDb = org.Mm.eg.db    # Mouse
```

### CellChat (Step 07)
```r
min.cells = 10          # L-R pair 발현 최소 세포 수
database = CellChatDB.human
signaling = "WNT"       # Pathway focus
```

---

## 📚 참고 방법론

### 핵심 알고리즘
- **Seurat:** Stuart et al. (2019). *Cell*
- **Harmony:** Korsunsky et al. (2019). *Nature Methods*
- **clusterProfiler:** Wu et al. (2021). *The Innovation*
- **CellChat:** Jin et al. (2021). *Nature Communications*

### 통계적 방법
- **DEG test:** Wilcoxon rank-sum (non-parametric)
- **Multiple testing correction:** Bonferroni
- **GO enrichment:** Hypergeometric test + BH adjustment

---

## 🎯 프로젝트 하이라이트

### 1. Multi-Omics Integration 능력
- Mouse와 Human 데이터 모두 다룸
- Species-specific 처리 방법 이해

### 2. 포괄적 분석 파이프라인
- QC부터 functional enrichment까지 end-to-end
- Cell-cell communication 분석 추가

### 3. 생물학적 통찰
- Immune microenvironment 특성화
- PORCN 결핍의 면역학적 효과 규명
- Cancer-immune crosstalk 메커니즘

### 4. 코드 품질
- Clean, readable, maintainable
- 면접 코드리뷰에 적합한 수준
- Production-ready structure

---

## 📝 기술적 Highlights

### Seurat v5 호환성
```r
# Layer management
seu <- JoinLayers(seu)  # Before FindAllMarkers

# Expression data access
expr <- GetAssayData(seu, slot = "data")
```

### Species-Specific 처리
```r
# Mouse
percent.mt <- PercentageFeatureSet(seu, pattern = "^mt-")

# Human (for CellChat)
rownames(data) <- toupper(str_remove(rownames(data), "\\.[0-9]+$"))
```

### Memory 최적화
```r
# Clean up
rm(unused_objects)
gc()
```

---

## 🤝 연락처

**분석 파이프라인 및 협업 문의:**
- Institution: Seoul National University Hospital, Mokam Research Institute
- Position: Computational Biologist (3.5 years experience)
- Specialization: Multi-omics analysis, scRNA-seq, CyTOF, IMC

---

## 📄 라이선스

이 분석 파이프라인은 교육 및 연구 목적으로 제공됩니다. 

---

## ⚠️ 중요 사항

### 데이터 구분
- **Steps 01-06:** Mouse PORCN KO vs WT (하나의 프로젝트)
- **Step 07:** Human Pancreatic Cancer (별도 프로젝트)

두 프로젝트는 독립적이며, 다양한 scRNA-seq 분석 역량을 보여주기 위해 함께 제시되었습니다.

### Raw Data
레포지토리에는 분석 코드만 포함되어 있습니다. Raw sequencing data는 크기 및 개인정보 보호 문제로 포함되지 않았습니다. 재현성을 위한 processed object는 요청시 제공 가능합니다.

---

**Last Updated:** 2025
**Pipeline Version:** 1.4
