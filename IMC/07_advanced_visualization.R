################################################################################
# IMC Analysis Pipeline
# Step 07: Advanced Visualization for Publication
# 
# Description:
#   - Effect size analysis (Cohen's d)
#   - Multi-panel comprehensive figure (Figure 1):
#     * Panel A: Category-wise violin + boxplot with statistics
#     * Panel B: Region-specific patterns (faceted)
#     * Panel C: Volcano plot (effect size vs p-value)
#     * Panel D: Top 15 distance changes by magnitude
#   - Distance change heatmap (Figure 2)
#   - Detailed comparison for key pairs (Figure 3)
#   - Supplementary tables with complete statistics
#
# Input:
#   - ./figure/06_spatial_analysis/region_distances_raw.csv
#   - ./figure/06_spatial_analysis/region_distances_statistics.csv
#
# Output:
#   - Publication-ready figures (TIFF 400-600 DPI + PDF)
#   - Supplementary data tables (CSV)
#
# Author: YMS
# Date: 2024-11-17
################################################################################

# Load required packages
library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggsci)
library(ggrepel)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Set parameters
ADV_VIZ_PARAMS <- list(
    volcano_label_threshold_d = 0.8,    # Label points with |d| > threshold
    volcano_label_threshold_p = 0.01,   # or p < threshold
    top_changes_n            = 15,      # Top N changes for panel D
    heatmap_top_n            = 30,      # Top N for heatmap
    key_pairs_per_category   = 2,       # Key pairs per category for Figure 3
    seed                     = 900223
)

# Create output directory
if (!dir.exists("./figure/07_publication")) {
    dir.create("./figure/07_publication", recursive = TRUE)
}

################################################################################
# 1. Load Distance Analysis Results
################################################################################

raw_data <- read_csv("./figure/06_spatial_analysis/region_distances_raw.csv")
stats_data <- read_csv("./figure/06_spatial_analysis/region_distances_statistics.csv")

cat("Data Overview:\n")
cat(sprintf("  Total measurements: %d\n", nrow(raw_data)))
cat(sprintf("  Combinations tested: %d\n", nrow(stats_data)))
cat(sprintf("  Samples: %d\n", length(unique(raw_data$sample_id))))
cat(sprintf("  Regions: %d\n", length(unique(raw_data$region))))
cat(sprintf("  Categories: %d\n\n", length(unique(raw_data$category))))

################################################################################
# 2. Effect Size Analysis (Cohen's d)
################################################################################

# IMPORTANT: Effect sizes are often more informative than p-values
# Cohen's d interpretation: 0.2=small, 0.5=medium, 0.8=large

effect_size_data <- raw_data %>%
    group_by(region, from_celltype, to_celltype, category, group) %>%
    summarize(
        mean_dist = mean(mean_dist, na.rm = TRUE),
        sd_dist = sd(mean_dist, na.rm = TRUE),
        n = n(),
        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from = group,
        values_from = c(mean_dist, sd_dist, n),
        names_sep = "_"
    ) %>%
    mutate(
        # Pooled standard deviation
        pooled_sd = sqrt(((n_Wt - 1) * sd_dist_Wt^2 + 
                          (n_Ko - 1) * sd_dist_Ko^2) / 
                         (n_Wt + n_Ko - 2)),
        
        # Cohen's d
        cohens_d = (mean_dist_Ko - mean_dist_Wt) / pooled_sd,
        
        # Absolute and relative change
        abs_change = mean_dist_Ko - mean_dist_Wt,
        pct_change = 100 * abs_change / mean_dist_Wt,
        
        # Effect size category
        effect_category = case_when(
            abs(cohens_d) < 0.2 ~ "Negligible",
            abs(cohens_d) < 0.5 ~ "Small",
            abs(cohens_d) < 0.8 ~ "Medium",
            TRUE ~ "Large"
        ),
        
        # Direction
        direction = ifelse(abs_change > 0, "Increased in Ko", "Decreased in Ko")
    ) %>%
    left_join(stats_data, by = c("region", "from_celltype", "to_celltype", "category"))

# Save effect size data
write_csv(
    effect_size_data,
    "./figure/07_publication/effect_size_analysis.csv"
)

cat("Effect Size Summary:\n")
cat(sprintf("  Large effects (|d| > 0.8): %d\n", 
            sum(abs(effect_size_data$cohens_d) > 0.8, na.rm = TRUE)))
cat(sprintf("  Medium effects (0.5 < |d| ≤ 0.8): %d\n",
            sum(abs(effect_size_data$cohens_d) > 0.5 & 
                abs(effect_size_data$cohens_d) <= 0.8, na.rm = TRUE)))
cat(sprintf("  Small effects (0.2 < |d| ≤ 0.5): %d\n\n",
            sum(abs(effect_size_data$cohens_d) > 0.2 & 
                abs(effect_size_data$cohens_d) <= 0.5, na.rm = TRUE)))

################################################################################
# 3. FIGURE 1: Multi-Panel Comprehensive Overview
################################################################################

cat("Creating Figure 1: Comprehensive overview...\n")

## Panel A: Category-wise Distance Comparison
category_data <- raw_data %>%
    group_by(category, sample_id, group) %>%
    summarize(mean_dist = mean(mean_dist, na.rm = TRUE), .groups = "drop")

category_stats <- category_data %>%
    group_by(category) %>%
    wilcox_test(mean_dist ~ group) %>%
    add_xy_position(x = "category")

p_a <- ggplot(category_data, aes(x = category, y = mean_dist, fill = group)) +
    geom_violin(alpha = 0.3, trim = FALSE) +
    geom_boxplot(
        width = 0.3,
        alpha = 0.7,
        position = position_dodge(0.9),
        outlier.shape = NA
    ) +
    geom_point(
        position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
        size = 1.5,
        alpha = 0.6
    ) +
    stat_pvalue_manual(
        category_stats,
        label = "p = {p}",
        tip.length = 0.01,
        size = 3
    ) +
    scale_fill_manual(
        values = c("Wt" = "#4DAF4A", "Ko" = "#E41A1C"),
        name = "Group"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_classic(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "black")
    ) +
    labs(
        x = "",
        y = "Mean Distance (μm)",
        title = "Cell-Cell Distance by Interaction Category"
    )

## Panel B: Region-wise Comparison
region_summary <- raw_data %>%
    group_by(region, category, group) %>%
    summarize(
        mean = mean(mean_dist, na.rm = TRUE),
        se = sd(mean_dist, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    )

p_b <- ggplot(region_summary, aes(x = region, y = mean, fill = group)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8, color = "black", size = 0.3) +
    geom_errorbar(
        aes(ymin = mean - se, ymax = mean + se),
        position = position_dodge(0.9),
        width = 0.3,
        size = 0.5
    ) +
    facet_wrap(~ category, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Wt" = "#4DAF4A", "Ko" = "#E41A1C")) +
    theme_classic(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom"
    ) +
    labs(
        x = "Region",
        y = "Mean Distance (μm)",
        fill = "Group",
        title = "Distance Patterns Across Regions"
    )

## Panel C: Volcano Plot (Effect Size vs P-value)
volcano_data <- effect_size_data %>%
    mutate(
        neg_log10_p = -log10(p),
        significant = ifelse(p < 0.05, "Nominal (p<0.05)", "NS"),
        large_effect = ifelse(abs(cohens_d) > 0.5, "Large Effect", "Small Effect"),
        label_text = ifelse(
            abs(cohens_d) > ADV_VIZ_PARAMS$volcano_label_threshold_d | 
            p < ADV_VIZ_PARAMS$volcano_label_threshold_p,
            paste0(from_celltype, "→", to_celltype, "\n(", region, ")"),
            ""
        )
    )

p_c <- ggplot(volcano_data, aes(x = cohens_d, y = neg_log10_p)) +
    annotate(
        "rect",
        xmin = -Inf, xmax = -0.5,
        ymin = -log10(0.05), ymax = Inf,
        fill = "#2166AC",
        alpha = 0.1
    ) +
    annotate(
        "rect",
        xmin = 0.5, xmax = Inf,
        ymin = -log10(0.05), ymax = Inf,
        fill = "#B2182B",
        alpha = 0.1
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        color = "grey40",
        size = 0.5
    ) +
    geom_vline(
        xintercept = c(-0.5, 0.5),
        linetype = "dashed",
        color = "grey40",
        size = 0.5
    ) +
    geom_point(aes(color = category, size = abs(abs_change)), alpha = 0.7) +
    geom_text_repel(
        aes(label = label_text),
        size = 2.5,
        max.overlaps = 15,
        box.padding = 0.5,
        segment.size = 0.3
    ) +
    scale_color_npg(name = "Category") +
    scale_size_continuous(name = "Absolute\nChange (μm)", range = c(1, 6)) +
    theme_classic(base_size = 11) +
    theme(legend.position = "right") +
    labs(
        x = "Effect Size (Cohen's d)\n← Decreased in Ko | Increased in Ko →",
        y = "-log10(p-value)",
        title = "Effect Size vs Statistical Significance"
    ) +
    annotate(
        "text",
        x = -1.5,
        y = max(volcano_data$neg_log10_p, na.rm = TRUE) * 0.95,
        label = "Decreased\n& Significant",
        size = 3,
        fontface = "italic",
        color = "#2166AC"
    ) +
    annotate(
        "text",
        x = 1.5,
        y = max(volcano_data$neg_log10_p, na.rm = TRUE) * 0.95,
        label = "Increased\n& Significant",
        size = 3,
        fontface = "italic",
        color = "#B2182B"
    )

## Panel D: Top Changes by Absolute Distance Change
top_changes_data <- effect_size_data %>%
    arrange(desc(abs(abs_change))) %>%
    head(ADV_VIZ_PARAMS$top_changes_n) %>%
    mutate(
        pair_label = paste0(from_celltype, " → ", to_celltype),
        pair_label = factor(pair_label, levels = pair_label)
    )

plot_data_d <- raw_data %>%
    semi_join(
        top_changes_data,
        by = c("region", "from_celltype", "to_celltype")
    ) %>%
    left_join(
        top_changes_data %>%
            select(region, from_celltype, to_celltype, pair_label),
        by = c("region", "from_celltype", "to_celltype")
    )

p_d <- ggplot(plot_data_d, aes(x = pair_label, y = mean_dist, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_point(
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
        size = 1,
        alpha = 0.5
    ) +
    scale_fill_manual(values = c("Wt" = "#4DAF4A", "Ko" = "#E41A1C")) +
    facet_wrap(~ region, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 9) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom"
    ) +
    labs(
        x = "",
        y = "Distance (μm)",
        fill = "Group",
        title = "Top 15 Distance Changes (by Magnitude)"
    )

## Combine panels
layout <- "
AABBBB
AABBBB
CCCCDD
CCCCDD
"

fig1 <- p_a + p_b + p_c + p_d +
    plot_layout(design = layout) +
    plot_annotation(
        title = "Spatial Cell-Cell Distance Analysis: Wt vs Ko",
        subtitle = "Region-based analysis of immune-myeloid-epithelial interactions",
        tag_levels = 'A',
        theme = theme(
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 12, color = "grey30")
        )
    )

ggsave(
    filename = "./figure/07_publication/FIGURE1_comprehensive.tiff",
    plot = fig1,
    width = 18,
    height = 14,
    dpi = 400,
    compression = "lzw"
)

ggsave(
    filename = "./figure/07_publication/FIGURE1_comprehensive.pdf",
    plot = fig1,
    width = 18,
    height = 14
)

cat("  ✓ Saved: FIGURE1_comprehensive (TIFF & PDF)\n\n")

################################################################################
# 4. FIGURE 2: Distance Change Heatmap
################################################################################

cat("Creating Figure 2: Distance change heatmap...\n")

# Select top N by effect size or p-value
heatmap_data <- effect_size_data %>%
    filter(!is.na(abs_change) & !is.na(cohens_d)) %>%
    arrange(desc(abs(cohens_d))) %>%
    head(ADV_VIZ_PARAMS$heatmap_top_n) %>%
    mutate(
        pair = paste0(from_celltype, " → ", to_celltype),
        region_pair = paste(region, pair, sep = ": ")
    )

# Create matrix
change_matrix <- heatmap_data %>%
    select(region, pair, abs_change) %>%
    pivot_wider(names_from = pair, values_from = abs_change) %>%
    column_to_rownames("region") %>%
    as.matrix()

# Column annotation (category)
category_annot <- heatmap_data %>%
    select(pair, category) %>%
    distinct() %>%
    arrange(pair)

col_ha <- HeatmapAnnotation(
    Category = category_annot$category,
    col = list(
        Category = c(
            "Immune->Epithelial" = "#E41A1C",
            "Immune->Myeloid" = "#377EB8",
            "Myeloid->Immune" = "#4DAF4A",
            "Myeloid->Epithelial" = "#FF7F00"
        )
    ),
    show_legend = TRUE,
    annotation_name_side = "left"
)

# Create heatmap
ht <- Heatmap(
    change_matrix,
    name = "Distance Change\n(Ko - Wt, μm)",
    
    col = colorRamp2(
        seq(
            -max(abs(change_matrix), na.rm = TRUE),
            max(abs(change_matrix), na.rm = TRUE),
            length.out = 9
        ),
        rev(brewer.pal(9, "RdBu"))
    ),
    
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 90,
    
    top_annotation = col_ha,
    
    column_title = sprintf("Top %d Cell-Cell Distance Changes (by Effect Size)",
                          ADV_VIZ_PARAMS$heatmap_top_n),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    
    na_col = "grey95",
    border = TRUE,
    
    heatmap_legend_param = list(
        legend_height = unit(4, "cm"),
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9)
    )
)

tiff(
    filename = "./figure/07_publication/FIGURE2_heatmap.tiff",
    width = 14,
    height = 9,
    units = "in",
    res = 400,
    compression = "lzw"
)
draw(ht)
dev.off()

pdf(
    filename = "./figure/07_publication/FIGURE2_heatmap.pdf",
    width = 14,
    height = 9
)
draw(ht)
dev.off()

cat("  ✓ Saved: FIGURE2_heatmap (TIFF & PDF)\n\n")

################################################################################
# 5. FIGURE 3: Detailed Comparison for Key Pairs
################################################################################

cat("Creating Figure 3: Detailed comparison for key pairs...\n")

# Select key pairs: largest effect sizes per category
key_pairs <- effect_size_data %>%
    group_by(category) %>%
    slice_max(abs(cohens_d), n = ADV_VIZ_PARAMS$key_pairs_per_category) %>%
    ungroup()

plot_list_f3 <- list()

for (i in 1:min(8, nrow(key_pairs))) {
    data_sub <- raw_data %>%
        filter(
            region == key_pairs$region[i],
            from_celltype == key_pairs$from_celltype[i],
            to_celltype == key_pairs$to_celltype[i]
        )
    
    if (nrow(data_sub) == 0) next
    
    p <- ggplot(data_sub, aes(x = group, y = mean_dist, fill = group)) +
        geom_violin(alpha = 0.3, trim = FALSE) +
        geom_boxplot(width = 0.3, alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
        scale_fill_manual(values = c("Wt" = "#4DAF4A", "Ko" = "#E41A1C")) +
        theme_classic(base_size = 10) +
        labs(
            title = sprintf(
                "%s\n%s → %s",
                key_pairs$region[i],
                key_pairs$from_celltype[i],
                key_pairs$to_celltype[i]
            ),
            subtitle = sprintf(
                "Effect size: %.2f | p = %.3f",
                key_pairs$cohens_d[i],
                key_pairs$p[i]
            ),
            x = "",
            y = "Distance (μm)"
        ) +
        stat_compare_means(
            method = "wilcox.test",
            label = "p.format",
            size = 3.5
        ) +
        theme(
            legend.position = "none",
            plot.title = element_text(size = 9, face = "bold"),
            plot.subtitle = element_text(size = 8, color = "grey30")
        )
    
    plot_list_f3[[i]] <- p
}

fig3 <- wrap_plots(plot_list_f3, ncol = 4) +
    plot_annotation(
        title = "Detailed Analysis: Cell Pairs with Largest Effect Sizes",
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

ggsave(
    filename = "./figure/07_publication/FIGURE3_key_pairs.tiff",
    plot = fig3,
    width = 16,
    height = 10,
    dpi = 400,
    compression = "lzw"
)

ggsave(
    filename = "./figure/07_publication/FIGURE3_key_pairs.pdf",
    plot = fig3,
    width = 16,
    height = 10
)

cat("  ✓ Saved: FIGURE3_key_pairs (TIFF & PDF)\n\n")

################################################################################
# 6. Supplementary Table: Complete Statistics
################################################################################

cat("Creating supplementary table...\n")

supp_table <- effect_size_data %>%
    select(
        Region = region,
        Category = category,
        From = from_celltype,
        To = to_celltype,
        Wt_mean = mean_dist_Wt,
        Ko_mean = mean_dist_Ko,
        Absolute_change = abs_change,
        Percent_change = pct_change,
        Cohens_d = cohens_d,
        Effect_category = effect_category,
        p_value = p,
        p_adjusted = p.adj,
        n_Wt,
        n_Ko
    ) %>%
    arrange(p_value) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)))

write_csv(
    supp_table,
    "./figure/07_publication/SUPPLEMENTARY_TABLE_complete_statistics.csv"
)

cat("  ✓ Saved: SUPPLEMENTARY_TABLE_complete_statistics.csv\n\n")

################################################################################
# 7. Summary Report
################################################################################

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("         PUBLICATION FIGURES COMPLETE                      \n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("Output Directory: ./figure/07_publication/\n\n")

cat("Generated Files:\n")
cat("───────────────────────────────────────────────────────────\n\n")

cat("Main Figures:\n")
cat("  • FIGURE1_comprehensive (TIFF/PDF)\n")
cat("    ├─ Panel A: Category comparison with statistics\n")
cat("    ├─ Panel B: Region-wise patterns\n")
cat("    ├─ Panel C: Volcano plot (effect size vs p-value)\n")
cat("    └─ Panel D: Top 15 distance changes\n\n")

cat("  • FIGURE2_heatmap (TIFF/PDF)\n")
cat("    └─ Top 30 distance changes by effect size\n\n")

cat("  • FIGURE3_key_pairs (TIFF/PDF)\n")
cat("    └─ Detailed analysis of key pairs\n\n")

cat("Data Tables:\n")
cat("  • effect_size_analysis.csv\n")
cat("    └─ Complete effect size calculations\n\n")

cat("  • SUPPLEMENTARY_TABLE_complete_statistics.csv\n")
cat("    └─ All statistical results with effect sizes\n\n")

cat("───────────────────────────────────────────────────────────\n\n")

# Summary statistics
cat("Key Findings:\n")
cat("───────────────────────────────────────────────────────────\n\n")

cat(sprintf("1. Total combinations analyzed: %d\n", nrow(effect_size_data)))

sig_count <- sum(effect_size_data$p.adj < 0.05, na.rm = TRUE)
if (sig_count > 0) {
    cat(sprintf("2. Significant pairs (p.adj < 0.05): %d\n", sig_count))
} else {
    cat("2. No significant differences after multiple testing correction\n")
    nominal_sig <- sum(effect_size_data$p < 0.05, na.rm = TRUE)
    cat(sprintf("   However, %d showed nominal significance (p<0.05)\n", nominal_sig))
}

large_effects <- sum(abs(effect_size_data$cohens_d) > 0.8, na.rm = TRUE)
cat(sprintf("\n3. Large effect sizes (|d| > 0.8): %d pairs\n", large_effects))

cat("\n───────────────────────────────────────────────────────────\n\n")

cat("Analysis complete!\n")
cat("All publication-ready figures and tables have been generated.\n\n")

# Clean up
rm(effect_size_data, volcano_data, heatmap_data, ht)
gc()

# NOTE: Pipeline complete!
