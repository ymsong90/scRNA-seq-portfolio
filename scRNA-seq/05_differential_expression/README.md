# Step 05: Differential Expression Analysis

## Overview
This module performs comprehensive differential gene expression (DEG) analysis comparing KO vs WT conditions across different cell populations. Includes cluster-specific comparisons, volcano plot generation, and statistical testing.

## Input Data
- `porcn.combined.harmony_annotated.RData`: Annotated Seurat object from Step 03

## Analysis Modules

### 05a: Myeloid Cell DEG Analysis (`05a_myeloid_DEG_volcano.R`)

**Purpose**: Identify differentially expressed genes in myeloid cell populations

**Key Features**:
- Cluster-specific KO vs WT comparisons
- Separate analysis for each myeloid subtype:
  - Classical Monocyte
  - Activated Monocyte
  - Trem2+ Macrophage
  - Selenop+ Macrophage
  - Hexb+ Macrophage
- Both raw and filtered DEG outputs
- Volcano plot visualization

**Parameters**:
```r
PARAMS <- list(
    deg = list(
        logfc = 0.25,        # Log fold-change threshold
        minpct = 0.10,       # Minimum detection rate
        padj = 0.05          # Adjusted p-value cutoff
    )
)
```

**Output**:
- `Myeloid_DEG_ko_vs_wt.csv`: All myeloid DEGs
- `DEG/`: Individual cluster DEG files
- `plots/volcano/`: Volcano plots for each cluster

### 05b: CD8 T Cell DEG Analysis (`05b_CD8T_DEG_analysis.R`)

**Purpose**: Analyze CD8 T cell transcriptional changes

**Two Comparisons**:
1. **CD8 Internal**: KO vs WT within CD8 T cells
2. **CD8 vs Rest**: CD8 T cells vs all other cell types

**Key Features**:
- Dual comparison strategy
- GO enrichment for up/down-regulated genes
- Custom volcano plot parameters per analysis
- DEG signature scoring

**Output**:
- `DEG_all/`: Raw DEG tables
- `DEG_up/`: Upregulated genes (padj < 0.05, logFC > 0.25)
- `DEG_down/`: Downregulated genes (padj < 0.05, logFC < -0.25)
- `GO_up/GO_down/`: Functional enrichment results
- `volcano/`: Customized volcano plots

## Volcano Plot Features

### Standard Elements
- **X-axis**: log2 fold change (KO vs WT)
- **Y-axis**: -log10(adjusted p-value)
- **Color coding**:
  - Red/Pink: Significantly upregulated
  - Blue/Purple: Significantly downregulated
  - Grey: Not significant

### Customization Options
```r
plot_volcano <- function(df, tag,
                        lfc_cut  = 0.25,    # Fold-change cutoff
                        padj_cut = 0.05,    # P-value cutoff
                        label_n  = 20)      # Number of genes to label
```

### Gene Labeling Strategy
- Top N most significant genes per direction (up/down)
- Uses `ggrepel` to avoid overlapping labels
- Adjustable via `label_n` parameter

## Statistical Testing

### Method
Wilcoxon rank-sum test (default in Seurat's `FindMarkers`)

### Why Wilcoxon?
- Non-parametric: No normality assumption
- Robust to outliers
- Standard in scRNA-seq field
- Good performance with sparse data

### Multiple Testing Correction
Bonferroni correction (Seurat default)

## DEG Interpretation Guidelines

### Significant DEG Criteria
1. |log2FC| > 0.25 (≈1.2-fold change)
2. p.adj < 0.05
3. Detected in ≥10% of cells in at least one group

### Biological Significance
- **log2FC > 1**: Strong upregulation (2-fold+)
- **log2FC 0.5-1**: Moderate upregulation
- **log2FC 0.25-0.5**: Mild upregulation

### Common Pitfalls
- **High logFC, low expression**: May not be biologically relevant
- **Low logFC, high significance**: Large sample size effect
- **Detection rate < 10%**: Unreliable estimates

## Usage

### Run Complete Myeloid Analysis
```r
source("05_differential_expression/05a_myeloid_DEG_volcano.R")
```

### Run CD8 T Cell Analysis
```r
source("05_differential_expression/05b_CD8T_DEG_analysis.R")
```

### Custom DEG Analysis
```r
# Load data
load("./data/porcn.combined.harmony_annotated.RData")

# Set identity to cell type of interest
Idents(porcn.combined.harmony) <- "NH_labels"

# Create condition × cell type metadata
porcn.combined.harmony$celltype.condition <- paste(
    Idents(porcn.combined.harmony),
    porcn.combined.harmony$ID,
    sep = "_"
)

# Run FindMarkers for specific comparison
Idents(porcn.combined.harmony) <- "celltype.condition"
my_deg <- FindMarkers(
    porcn.combined.harmony,
    ident.1 = "YourCellType_KO",
    ident.2 = "YourCellType_WT",
    min.pct = 0.1,
    logfc.threshold = 0.25
)
```

## Output File Structure

```
05_differential_expression/
├── results/
│   ├── DEG/
│   │   ├── Myeloid_DEG_ko_vs_wt.csv
│   │   ├── CD8_KO_vs_WT.csv
│   │   └── CD8_vs_Rest.csv
│   ├── DEG_up/
│   │   ├── DEG_Up_CD8_KO_vs_WT.csv
│   │   └── ...
│   ├── DEG_down/
│   │   └── ...
│   └── plots/
│       └── volcano/
│           ├── Volcano_Classical_Monocyte.tiff
│           ├── Volcano_CD8_KO_vs_WT.tiff
│           └── ...
```

## Expected Results

### Typical DEG Counts per Cluster
- **Monocytes**: 500-2,000 DEGs
- **Macrophages**: 1,000-3,000 DEGs
- **T cells**: 300-1,500 DEGs

### Top Gene Categories Often Affected
1. **Immune response genes**
2. **Metabolic genes**
3. **Cell cycle genes** (in proliferative populations)
4. **Chemokines/cytokines**

## Troubleshooting

### Issue: No significant DEGs found
**Solutions**:
1. Lower logfc.threshold to 0.1
2. Increase min.pct to capture more genes
3. Check if groups have sufficient cells (≥50 recommended)

### Issue: Too many DEGs (>5,000)
**Solutions**:
1. Increase stringency: logfc.threshold = 0.5, padj = 0.01
2. Check for batch effects
3. Verify cell type annotations

### Issue: Volcano plots are crowded
**Solutions**:
1. Reduce `label_n` parameter
2. Increase `max.overlaps` in `geom_text_repel`
3. Filter genes before plotting (e.g., only high-confidence genes)

## Performance Optimization

### For Large Datasets
```r
# Use parallel processing
library(future)
plan(multisession, workers = 8)

# Use future.apply for cluster loops
cluster_degs <- future_lapply(cluster_ids, function(cl) {
    FindMarkers(...)
}, future.seed = TRUE)
```

### Memory Management
```r
# Clear unnecessary objects
rm(intermediate_objects)
gc()

# Use logfc.threshold = 0 only when necessary
# Pre-filtering saves memory and time
```

## Dependencies

- Seurat (≥5.0)
- dplyr, tibble
- ggplot2, ggrepel
- future, future.apply (for parallel processing)

## Integration with GO Analysis

DEG results automatically flow into Step 06 (Functional Enrichment):
```r
# Example workflow
source("05_differential_expression/05a_myeloid_DEG_volcano.R")
source("06_functional_enrichment/06a_myeloid_GO_analysis.R")
```

## Next Step
Functional enrichment analysis: `06_functional_enrichment/`
