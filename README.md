# Single-Cell RNA-seq Analysis Pipeline

Comprehensive scRNA-seq analysis pipeline for PORCN KO vs WT comparison study, focusing on immune cell profiling and differential expression analysis.

## 📊 Project Overview

This repository contains a complete single-cell RNA sequencing analysis pipeline used to characterize the immune microenvironment in PORCN knockout vs wild-type conditions. The analysis identifies cell type-specific transcriptional changes and functional enrichments across multiple immune populations.

**Key Findings:**
- Comprehensive profiling of 20+ distinct cell populations
- Myeloid cell subset characterization with functional annotations
- Condition-specific differential gene expression across major immune lineages
- Pathway enrichment analysis revealing biological processes altered by PORCN deficiency

## 🔬 Technical Stack

- **Language:** R (version ≥ 4.0)
- **Core Packages:**
  - `Seurat` (v5+): Single-cell analysis framework
  - `Harmony`: Batch effect correction
  - `clusterProfiler`: GO/KEGG enrichment analysis
  - `ggplot2`, `ComplexHeatmap`: Visualization
- **Analysis Scope:**
  - ~20,000-50,000 cells analyzed
  - 2 conditions (WT vs KO)
  - Multiple cell type annotations and subclustering

## 📁 Repository Structure

```
scRNA-seq-portfolio/
├── README.md                          # This file
├── 01_preprocessing/                  # Data loading and QC
│   ├── 01_data_loading_QC.R
│   └── README.md
├── 02_integration/                    # Batch correction
│   ├── 02_harmony_integration.R
│   └── README.md
├── 03_clustering/                     # Clustering and annotation
│   ├── 03_clustering_annotation.R
│   ├── 03b_annotation_revision.R
│   └── README.md
├── 04_visualization/                  # Cell proportion analysis
│   ├── 04_cell_proportion_analysis.R
│   └── README.md
├── 05_differential_expression/        # DEG analysis
│   ├── 05a_myeloid_DEG_analysis.R
│   ├── 05b_CD8T_DEG_analysis.R
│   ├── 05c_volcano_plots.R
│   └── README.md
├── 06_functional_enrichment/          # GO/KEGG analysis
│   ├── 06a_myeloid_GO_analysis.R
│   ├── 06b_monocyte_macrophage_GO.R
│   ├── 06c_wnt_pathway_analysis.R
│   └── README.md
└── utils/                             # Helper functions
    ├── plotting_functions.R
    ├── deg_functions.R
    └── go_functions.R
```

## 🚀 Analysis Pipeline

### 1. Data Preprocessing & QC
```r
# Load 10X data
# Quality control (nFeature, nCount, percent.mt)
# Filter cells: 200 < nFeature < 8000, percent.mt < 20
```
**Output:** Filtered Seurat object

### 2. Normalization & Integration
```r
# LogNormalize (scale factor: 10,000)
# Find variable features (top 2,000)
# Harmony batch correction between WT/KO
```
**Output:** Batch-corrected integrated object

### 3. Dimensionality Reduction & Clustering
```r
# PCA (top 10 PCs)
# UMAP/t-SNE visualization
# Graph-based clustering (Louvain)
# Cell type annotation using canonical markers
```
**Output:** Annotated Seurat object with 20+ cell types

### 4. Cell Type Proportion Analysis
```r
# Calculate cell type distributions
# WT vs KO proportion comparisons
# Statistical testing for compositional changes
```
**Output:** Proportion plots and statistics

### 5. Differential Expression Analysis
```r
# Cluster-specific KO vs WT comparisons
# FindMarkers (Wilcoxon rank-sum test)
# Filtering: |log2FC| > 0.25, p.adj < 0.05
```
**Output:** DEG tables and volcano plots

### 6. Functional Enrichment
```r
# GO Biological Process enrichment
# KEGG pathway analysis
# Visualization: bar plots, dot plots
```
**Output:** Enriched terms and pathway visualizations

## 📈 Key Analyses

### Major Cell Populations Identified
- **Epithelial/Cancer cells** (6 subtypes)
- **T cells**: CD4 T, CD8 T (CTL), NK cells
- **Myeloid cells**: 
  - Monocytes (Classical, Activated)
  - Macrophages (Trem2+, Selenop+, Hexb+)
  - Dendritic cells (Batf3hi DC)
  - Neutrophils
- **B cells** (3 subtypes)
- **Stromal cells**: CAFs (4 subtypes), Endothelial cells

### Myeloid Subclustering Analysis
Focused analysis on myeloid compartment revealing:
- 5 distinct monocyte/macrophage populations
- Differential activation states between WT/KO
- Wnt pathway-related gene expression patterns

### Condition-Specific DEG Analysis
- Per-cluster KO vs WT comparisons
- CD8 T cell activation state analysis
- Pathway enrichment for up/down-regulated genes

## 💻 Usage

### Prerequisites
```r
# Install required packages
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"))
devtools::install_github("immunogenomics/harmony")
```

### Running the Pipeline
```r
# 1. Start with preprocessing
source("01_preprocessing/01_data_loading_QC.R")

# 2. Run integration
source("02_integration/02_harmony_integration.R")

# 3. Perform clustering
source("03_clustering/03_clustering_annotation.R")

# 4. Cell proportion analysis
source("04_visualization/04_cell_proportion_analysis.R")

# 5. Differential expression
source("05_differential_expression/05a_myeloid_DEG_analysis.R")

# 6. Functional enrichment
source("06_functional_enrichment/06a_myeloid_GO_analysis.R")
```

### Input Data Structure
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

## 📊 Example Outputs

### UMAP Visualization
Cell type annotations displayed on UMAP projection with condition (WT/KO) split view.

### Cell Proportion Stacked Bar Plots
Comparative analysis showing relative abundance of each cell type between conditions.

### Volcano Plots
Differential gene expression visualized for each cluster comparison (KO vs WT).

### GO Enrichment Bar/Dot Plots
Top enriched biological processes for upregulated and downregulated genes.

## 🔍 Key Features

### Code Quality
- **Modular structure**: Each analysis step in separate, well-documented scripts
- **Reproducibility**: Fixed random seeds, explicit parameters
- **Error handling**: Input validation and informative error messages
- **Performance**: Parallelization where applicable (future.apply)

### Visualization Standards
- Publication-quality figures (TIFF, 600 DPI)
- Consistent color schemes across plots
- Clear axis labels and legends
- Customizable themes

### Documentation
- Inline comments explaining complex logic
- Parameter documentation
- Expected input/output descriptions
- Troubleshooting guides in each module README

## 📝 Analysis Parameters

### Quality Control
- Minimum features per cell: 200
- Maximum features per cell: 8,000
- Maximum mitochondrial percentage: 20%

### Differential Expression
- Log fold-change threshold: 0.25
- Minimum detection rate: 10% (min.pct)
- Adjusted p-value cutoff: 0.05
- Statistical test: Wilcoxon rank-sum

### GO Enrichment
- Ontology: Biological Process (BP)
- P-value cutoff: 0.10
- Q-value cutoff: 0.20
- Organism: *Mus musculus*

## 🔬 Biological Insights

### PORCN Deficiency Effects
1. **Immune Cell Composition Changes**
   - Altered myeloid cell proportions
   - Shifts in T cell activation states

2. **Transcriptional Reprogramming**
   - Wnt pathway genes significantly affected
   - Immune response genes dysregulated

3. **Cell Type-Specific Responses**
   - Macrophages show strongest transcriptional changes
   - T cells exhibit altered effector profiles

## 📚 References

Key methodologies implemented:
- Stuart et al. (2019). Comprehensive Integration of Single-Cell Data. *Cell*
- Korsunsky et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*
- Wu et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*

## 🤝 Contact

For questions about the analysis pipeline or collaboration opportunities:
- GitHub: [Your GitHub Profile]
- Email: [Your Email]
- Institution: Seoul National University Hospital, Mokam Research Institute

## 📄 License

This analysis pipeline is provided for educational and research purposes. If you use this code, please cite appropriately.

---

**Note:** This repository contains analysis code only. Raw sequencing data is not included due to size and privacy considerations. Example processed objects may be provided upon request for reproducibility purposes.
