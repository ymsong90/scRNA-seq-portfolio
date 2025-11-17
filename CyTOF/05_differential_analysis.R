################################################################################
# Step 05: Differential Analysis
#
# Purpose: Differential abundance and state analysis between conditions
# Dataset: NMIBC bladder cancer CyTOF
#
# Input:
#   - ./data/sce_annotated.RData
#   - ./data/sce_myeloid_annotated.RData
#
# Output:
#   - ./results/05_differential/*.csv (statistical results)
#   - ./results/05_differential/*.png (visualizations)
#
# Author: YMS
# Date: 2025
################################################################################

library(CATALYST)
library(dplyr)
library(ggplot2)
library(tidyr)
library(diffcyt)

# Create output directory
if (!dir.exists("./results/05_differential")) {
    dir.create("./results/05_differential", recursive = TRUE)
}

################################################################################
# 1. Load Annotated Data
################################################################################

load("./data/sce_annotated.RData")

################################################################################
# 2. Differential Abundance (DA) Analysis
################################################################################

# NOTE: DA tests identify cell populations that change in frequency
# between conditions (e.g., BCG responders vs non-responders)

# Run differential abundance test
da_results <- diffcyt(
    sce_fil,
    experiment_info = metadata(sce_fil)$experiment_info,
    marker_info = rowData(sce_fil),
    design = ~ condition,
    contrast = c(0, 1, -1, 0, 0),  # Example: RE vs UR comparison
    analysis_type = "DA",
    method_DA = "diffcyt-DA-edgeR",
    min_cells = 3,
    min_samples = 2
)

# Extract and save results
da_table <- rowData(da_results$res)
write.csv(
    da_table,
    file = "./results/05_differential/DA_results_RE_vs_UR.csv",
    row.names = FALSE
)

# Print significant results
sig_da <- da_table[da_table$p_adj < 0.05, ]
print("Significant DA clusters:")
print(sig_da[, c("cluster_id", "logFC", "p_val", "p_adj")])

################################################################################
# 3. Differential State (DS) Analysis
################################################################################

# NOTE: DS tests identify markers with differential expression
# within specific cell types

ds_results <- diffcyt(
    sce_fil,
    experiment_info = metadata(sce_fil)$experiment_info,
    marker_info = rowData(sce_fil),
    design = ~ condition,
    contrast = c(0, 1, -1, 0, 0),
    analysis_type = "DS",
    method_DS = "diffcyt-DS-limma",
    markers_to_test = rowData(sce_fil)$marker_name[
        rowData(sce_fil)$marker_class == "state"
    ]
)

# Extract and save results
ds_table <- rowData(ds_results$res)
write.csv(
    ds_table,
    file = "./results/05_differential/DS_results_RE_vs_UR.csv",
    row.names = FALSE
)

# Print significant results
sig_ds <- ds_table[ds_table$p_adj < 0.05, ]
print("Significant DS markers:")
print(sig_ds[, c("cluster_id", "marker_id", "logFC", "p_val", "p_adj")])

################################################################################
# 4. Marker Expression Comparison
################################################################################

# NOTE: Compare median marker expression between conditions
# Useful for identifying functional state changes

# Define markers for comparison
comparison_markers <- c(
    'CD45', 'CD3e', 'CD4', 'CD8', 'CD11b', 'CD11c',
    'CD14', 'CD15', 'CD16', 'CD20', 'CD32', 'CD36',
    'CD38', 'CD64', 'CD68', 'CD192'
)

# Extract median expression data
pb_data <- plotPbExprs(
    sce_fil,
    k = "celltype",
    fun = "median",
    features = comparison_markers,
    facet_by = "cluster_id"
)

# Calculate fold changes
fold_change_data <- pb_data$data %>%
    group_by(cluster_id, antigen) %>%
    pivot_wider(names_from = condition, values_from = value) %>%
    mutate(
        fold_change = RE / UR,
        log2_fold_change = log2(RE / UR)
    )

write.csv(
    fold_change_data,
    file = "./results/05_differential/marker_fold_changes_RE_vs_UR.csv",
    row.names = FALSE
)

# Visualize fold changes
fc_plot <- ggplot(
    fold_change_data,
    aes(x = antigen, y = log2_fold_change, fill = log2_fold_change > 0)
) +
    geom_bar(stat = "identity", width = 0.7) +
    facet_wrap(~ cluster_id, scales = "free_x") +
    labs(
        x = "Marker",
        y = "Log2 Fold Change (RE / UR)",
        title = "Marker Expression Changes: RE vs UR"
    ) +
    scale_fill_manual(
        values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
        labels = c("Upregulated", "Downregulated"),
        name = "Direction"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30")

ggsave(
    filename = "./results/05_differential/marker_fold_changes_barplot.png",
    plot = fc_plot,
    width = 18,
    height = 12,
    dpi = 300
)

################################################################################
# 5. Line Plot: Marker Expression Across Conditions
################################################################################

# Extract data for line plots
line_data <- pb_data$data %>%
    group_by(cluster_id, antigen, condition) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

write.csv(
    line_data,
    file = "./results/05_differential/marker_expression_by_condition.csv",
    row.names = FALSE
)

# Generate line plot
line_plot <- ggplot(
    line_data,
    aes(x = antigen, y = mean_value, group = condition, color = condition)
) +
    geom_line(size = 0.8) +
    geom_point(size = 2) +
    labs(
        x = "Markers",
        y = "Mean Expression (arcsinh-transformed)",
        color = "Condition",
        title = "Marker Expression Across Conditions"
    ) +
    facet_wrap(~ cluster_id, scales = "free_y") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top",
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid = element_blank()
    )

ggsave(
    filename = "./results/05_differential/marker_expression_lineplot.png",
    plot = line_plot,
    width = 16,
    height = 12,
    dpi = 300
)

################################################################################
# 6. Myeloid Subtype-Specific Analysis
################################################################################

# NOTE: Analyze marker expression in myeloid cell subtypes
# Focus on activation and polarization markers

load("./data/sce_myeloid_annotated.RData")

# Define myeloid-specific markers
myeloid_markers <- c(
    'CD14', 'CD16', 'CD64', 'CD68', 'CD163', 'CD169',
    'CD206', 'HLA-DR', 'CD86', 'CD40', 'CD274',
    'Tim-3', 'VISTA', 'CD39', 'CD32'
)

# Extract myeloid marker expression
myeloid_pb <- plotPbExprs(
    sce_myeloid,
    k = "myeloid_subtype",
    fun = "median",
    features = myeloid_markers,
    facet_by = "cluster_id"
)

# Calculate fold changes for myeloid subtypes
myeloid_fc <- myeloid_pb$data %>%
    group_by(cluster_id, antigen) %>%
    pivot_wider(names_from = condition, values_from = value) %>%
    mutate(
        fold_change = RE / UR,
        log2_fold_change = log2(RE / UR)
    )

write.csv(
    myeloid_fc,
    file = "./results/05_differential/myeloid_fold_changes.csv",
    row.names = FALSE
)

# Visualize myeloid fold changes
myeloid_fc_plot <- ggplot(
    myeloid_fc,
    aes(x = antigen, y = log2_fold_change, fill = log2_fold_change > 0)
) +
    geom_bar(stat = "identity", width = 0.7) +
    facet_wrap(~ cluster_id, scales = "free_x") +
    labs(
        x = "Marker",
        y = "Log2 Fold Change (RE / UR)",
        title = "Myeloid Marker Expression Changes"
    ) +
    scale_fill_manual(
        values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
        name = "Direction"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = "dashed")

ggsave(
    filename = "./results/05_differential/myeloid_fold_changes_barplot.png",
    plot = myeloid_fc_plot,
    width = 14,
    height = 10,
    dpi = 300
)

################################################################################
# 7. Statistical Summary
################################################################################

# Create summary statistics
summary_stats <- data.frame(
    Analysis = c(
        "DA Test (Clusters)",
        "DA Significant (p_adj < 0.05)",
        "DS Test (Markers)",
        "DS Significant (p_adj < 0.05)"
    ),
    Count = c(
        nrow(da_table),
        sum(da_table$p_adj < 0.05, na.rm = TRUE),
        nrow(ds_table),
        sum(ds_table$p_adj < 0.05, na.rm = TRUE)
    )
)

write.csv(
    summary_stats,
    file = "./results/05_differential/analysis_summary.csv",
    row.names = FALSE
)

print(summary_stats)

################################################################################
# 8. Clean Up
################################################################################

rm(da_results, ds_results, pb_data, myeloid_pb)
gc()

################################################################################
# End of Step 05
################################################################################
