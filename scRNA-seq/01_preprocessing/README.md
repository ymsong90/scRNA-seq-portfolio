# Step 01: Data Loading and Quality Control

## Overview
This module handles the initial loading of 10X Genomics single-cell RNA-seq data and performs comprehensive quality control filtering to ensure high-quality downstream analysis.

## Input Data
- **Format:** 10X Genomics filtered feature-barcode matrices
- **Structure:**
  ```
  data/
  ├── wt/filtered_feature_bc_matrix/
  │   ├── barcodes.tsv.gz
  │   ├── features.tsv.gz
  │   └── matrix.mtx.gz
  └── ko/filtered_feature_bc_matrix/
      ├── barcodes.tsv.gz
      ├── features.tsv.gz
      └── matrix.mtx.gz
  ```

## Key Steps

### 1. Data Loading
- Reads 10X data using `Seurat::Read10X()`
- Creates separate Seurat objects for WT and KO conditions
- Initial filtering: genes detected in ≥3 cells, cells with ≥200 genes

### 2. Sample Merging
- Combines WT and KO datasets into a single Seurat object
- Adds sample ID metadata for tracking

### 3. QC Metric Calculation
- **nFeature_RNA**: Number of genes detected per cell
- **nCount_RNA**: Total UMI counts per cell
- **percent.mt**: Percentage of mitochondrial gene expression

### 4. Quality Filtering
Default thresholds (adjustable in script):
- Minimum features: 200 genes
- Maximum features: 8,000 genes
- Maximum mitochondrial percentage: 20%

Rationale:
- **Low nFeature**: Likely empty droplets or low-quality cells
- **High nFeature**: Possible doublets
- **High percent.mt**: Stressed/dying cells

### 5. Visualization
Generates before/after QC plots:
- Violin plots showing distribution of QC metrics
- Scatter plots examining relationships between metrics

## Output Files

### R Objects
- `porcn.combined_QC.RData`: QC-filtered Seurat object

### Figures
- `QC_violin_plots_before_filtering.png`: Initial QC metrics
- `QC_violin_plots_after_filtering.png`: Post-filtering QC metrics
- `QC_scatter_plots.png`: Metric relationships

## Usage

```r
# Run the complete QC pipeline
source("01_preprocessing/01_data_loading_QC.R")

# Or load the processed object directly
load("./data/porcn.combined_QC.RData")
```

## Expected Results

Typical filtering statistics:
- **Cell retention**: ~85-95% of input cells
- **Median nFeature**: 2,000-4,000 genes/cell
- **Median percent.mt**: 2-8%

## Troubleshooting

### Issue: Too many cells filtered out (>30%)
**Solution**: Check QC thresholds - may be too stringent
```r
QC_PARAMS <- list(
    min_features = 200,    # Try lowering
    max_features = 10000,  # Try increasing
    max_mt_pct   = 25      # Try increasing
)
```

### Issue: High mitochondrial percentage across all cells
**Solution**: May indicate sample quality issues - review experimental protocol

### Issue: Bimodal distribution in nFeature
**Solution**: May indicate presence of multiple cell populations or doublets

## Parameters

All parameters are defined at the top of the script:

```r
QC_PARAMS <- list(
    min_features = 200,      # Minimum genes per cell
    max_features = 8000,     # Maximum genes per cell
    max_mt_pct   = 20        # Maximum mitochondrial percentage
)
```

## Dependencies

- Seurat (≥5.0)
- dplyr
- ggplot2
- patchwork
- hdf5r

## Next Step
Proceed to normalization and integration: `02_integration/02_harmony_integration.R`
