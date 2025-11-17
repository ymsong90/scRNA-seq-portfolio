# 코딩 면접 준비 가이드

## 📌 면접 형식 대응 전략

### 1. GitHub를 통한 코드 베이스라인 리뷰

#### 준비 사항
✅ **Repository를 Public으로 설정**
- Settings → Change repository visibility → Make public

✅ **README.md 최적화**
- 프로젝트 개요를 명확하게
- 기술 스택 명시
- 주요 분석 결과 요약
- 시각적 요소 포함 (UMAP, Volcano plot 등)

✅ **코드 정리 체크리스트**
- [ ] 모든 파일에 주석이 충분한가?
- [ ] 변수명이 명확한가?
- [ ] 함수가 모듈화되어 있는가?
- [ ] 하드코딩된 경로가 없는가?
- [ ] 에러 처리가 적절한가?

#### Repository 링크 제공 방법
```
안녕하세요,

제 single-cell RNA-seq 분석 포트폴리오는 다음 GitHub 저장소에서 확인하실 수 있습니다:

https://github.com/[your-username]/scRNA-seq-portfolio

주요 분석 내용:
• 데이터 전처리 및 QC (01_preprocessing/)
• Harmony batch correction (02_integration/)
• Cell type annotation (03_clustering/)
• Differential expression analysis (05_differential_expression/)
• Functional enrichment (06_functional_enrichment/)

각 디렉토리의 README 파일에 상세한 분석 방법과 결과가 기술되어 있습니다.
실행 가능한 전체 코드와 함께 분석 파이프라인을 제공하고 있습니다.

감사합니다.
```

---

## 🎯 예상 질문 & 답변 준비

### 기술적 질문

#### Q1: "Harmony를 사용한 이유는 무엇인가요? CCA/RPCA와 비교했을 때 장점은?"
**답변**:
- Harmony는 linear correction 기반으로 빠르고 효율적
- WT vs KO 두 조건만 있어 Harmony가 충분
- Batch effect가 크지 않아 복잡한 방법 불필요
- 실제로 convergence plot으로 좋은 integration 확인

#### Q2: "QC 필터링 기준을 어떻게 정했나요?"
**답변**:
```r
# 근거 제시
- nFeature 200-8000: 
  * 200 미만: empty droplets
  * 8000 초과: potential doublets
- percent.mt < 20%: 
  * 높은 mt%는 stressed/dying cells
  * Mouse 데이터 기준 20%가 일반적
```

#### Q3: "FindMarkers에서 Wilcoxon test를 사용한 이유는?"
**답변**:
- Non-parametric test: 정규성 가정 불필요
- scRNA-seq 데이터는 zero-inflated하고 non-normal
- Field standard (Seurat default)
- Robust to outliers

#### Q4: "Multiple testing correction은 어떻게 했나요?"
**답변**:
- Bonferroni correction (Seurat default)
- Conservative하지만 false positive 최소화
- p.adj < 0.05로 엄격하게 필터링

### 분석 설계 질문

#### Q5: "Myeloid 세포만 따로 subset해서 분석한 이유는?"
**답변**:
- 전체 데이터에서는 major population만 구분
- Myeloid 내부 heterogeneity 파악 위해 re-clustering
- Monocyte/Macrophage subtype 구분 필요
- Finer resolution clustering 가능

#### Q6: "GO 분석에서 universe를 지정한 이유는?"
**답변**:
```r
# Background genes (universe) 지정의 중요성
bg_entrez <- bitr(rownames(seu), ...)  # 실제 측정된 유전자만
# 이유:
# 1. 정확한 통계: 실제 측정 가능한 유전자 대비 enrichment
# 2. Bias 제거: 측정되지 않은 유전자는 제외
# 3. False positive 감소
```

### 문제 해결 질문

#### Q7: "분석 중 어려웠던 점과 해결 방법은?"
**답변 예시**:
1. **Batch effect 문제**
   - 문제: 초기 UMAP에서 condition별로 분리됨
   - 해결: Harmony integration 적용
   - 검증: Integration plot으로 확인

2. **메모리 부족**
   - 문제: FindAllMarkers에서 메모리 에러
   - 해결: JoinLayers() 먼저 실행, gc() 활용
   - 최적화: future.apply로 병렬 처리

3. **Annotation 어려움**
   - 문제: 일부 cluster의 identity 불분명
   - 해결: Canonical marker 여러 개 조합 확인
   - 검증: Published literature와 비교

---

## 💻 실시간 코딩 면접 대비

### 준비할 코드 스니펫

#### 1. 빠른 QC 체크
```r
# 데이터 로딩
seu <- readRDS("your_data.rds")

# QC metrics 계산
seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^mt-")

# 분포 확인
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# Filtering
seu_filtered <- subset(seu, 
    subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
```

#### 2. Marker 찾기
```r
# 특정 cluster의 marker genes
Idents(seu) <- "seurat_clusters"
markers <- FindMarkers(
    seu,
    ident.1 = "5",  # Target cluster
    min.pct = 0.25,
    logfc.threshold = 0.25
)

# Top markers
head(markers %>% arrange(p_val_adj), 10)
```

#### 3. Visualization
```r
# Feature plot
FeaturePlot(seu, features = c("Cd68", "Cd3e", "Krt19"))

# DotPlot
DotPlot(seu, features = c("Cd68", "Adgre1", "Itgam"))

# UMAP with custom colors
DimPlot(seu, group.by = "celltype", cols = my_colors)
```

#### 4. DEG 분석
```r
# Condition 비교
seu$celltype.condition <- paste(seu$celltype, seu$condition, sep = "_")
Idents(seu) <- "celltype.condition"

deg <- FindMarkers(
    seu,
    ident.1 = "Macrophage_KO",
    ident.2 = "Macrophage_WT",
    logfc.threshold = 0.25,
    min.pct = 0.1
)

# Significant genes
sig_genes <- deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)
```

---

## 🎭 시연 시나리오

### 5분 Quick Demo
```r
# 1. 데이터 로딩 (30초)
load("porcn.combined.harmony_annotated.RData")

# 2. 전체 UMAP 보기 (30초)
DimPlot(porcn.combined.harmony, group.by = "NH_labels", label = TRUE)

# 3. 특정 cell type marker 확인 (1분)
FeaturePlot(porcn.combined.harmony, 
    features = c("Cd68", "Cd3e", "Krt19", "Pecam1"))

# 4. DEG 결과 확인 (2분)
myeloid_deg <- read.csv("results/DEG/Myeloid_DEG_ko_vs_wt.csv")
head(myeloid_deg %>% arrange(p_val_adj), 20)

# Volcano plot
source("utils/plotting_functions.R")
p <- plot_volcano_custom(myeloid_deg, title = "Myeloid: KO vs WT")
print(p)

# 5. GO 결과 설명 (1분)
go_results <- read.csv("results/GO/GO_Macrophage_up.csv")
head(go_results, 10)
```

### 15분 상세 Demo
위 내용 + 추가:
- Subset 분석 과정 설명
- 코드 구조 설명
- 파라미터 선택 이유 설명
- 결과 해석 및 생물학적 의미 논의

---

## 📋 체크리스트: 면접 전날

### GitHub Repository
- [ ] Public으로 설정됨
- [ ] README.md 최신화
- [ ] 모든 코드 파일에 주석 충분
- [ ] 불필요한 파일 제거 (.gitignore 확인)
- [ ] 예시 결과 figure 포함

### 로컬 환경
- [ ] 모든 스크립트 실행 테스트 완료
- [ ] 주요 결과 파일 확인
- [ ] Demo용 데이터 준비
- [ ] 필요한 패키지 모두 설치됨

### 설명 준비
- [ ] 분석 전체 플로우 설명 가능
- [ ] 각 단계별 선택 이유 설명 가능
- [ ] 주요 결과 해석 준비
- [ ] 어려웠던 점과 해결법 준비

### 예상 질문 답변
- [ ] 기술적 질문 답변 준비
- [ ] 생물학적 해석 준비
- [ ] Alternative approach 생각해봄
- [ ] Limitation & Future work 정리

---

## 💡 면접 팁

### Do's ✅
1. **코드 설명 시**
   - "이 부분은 ~를 위한 것입니다"
   - "~한 이유로 이 방법을 선택했습니다"
   - "결과는 ~를 의미합니다"

2. **질문 받을 때**
   - 이해 못하면 다시 물어보기
   - 생각할 시간 요청하기
   - 확신 없으면 "~라고 생각하는데, 확실하지 않습니다" 솔직히 말하기

3. **코드 작성 시**
   - 주석 먼저 작성 (의사코드)
   - 단계별로 테스트
   - 변수명 명확하게

### Don'ts ❌
1. 모르는 것을 아는 척
2. 너무 빨리 대답 (생각 없이)
3. 면접관 말 자르기
4. 방어적 태도

### 예상치 못한 질문 대처
```
"잠깐 생각할 시간을 주시겠습니까?"
"비슷한 경험은 ~입니다만, 정확히는 시도해보지 않았습니다"
"이론적으로는 ~하면 될 것 같은데, 검증은 필요할 것 같습니다"
```

---

## 🎓 추가 학습 리소스

### 복습 필수 개념
1. **통계**
   - Wilcoxon test vs t-test
   - Multiple testing correction
   - P-value vs FDR

2. **scRNA-seq**
   - Normalization methods
   - Batch correction
   - Clustering algorithms

3. **R Programming**
   - dplyr pipelines
   - ggplot2 syntax
   - Seurat functions

### 관련 논문 리뷰
- Seurat v5 paper
- Harmony paper
- clusterProfiler paper

---

**Good luck with your interview! 🍀**
