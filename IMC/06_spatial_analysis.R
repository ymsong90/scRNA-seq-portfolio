################################################################################
# IMC Analysis Pipeline
# Step 06: Spatial Distance Analysis
# 
# Description:
#   - Define cell type categories (Epithelial, Immune, Myeloid)
#   - Region-based cell-cell distance analysis
#   - Calculate minimum distances between cell type pairs
#   - Four interaction categories:
#     * Immune → Epithelial (immune infiltration to tumor)
#     * Immune → Myeloid (T cell proximity to myeloid)
#     * Myeloid → Immune (myeloid proximity to T cells)
#     * Myeloid → Epithelial (myeloid infiltration to tumor)
#   - Statistical testing (Wilcoxon rank-sum test: Wt vs Ko)
#   - Multiple testing correction (Benjamini-Hochberg FDR)
#
# Input:
#   - ./data/spe_filtered.rds (from Step 04)
#
# Output:
#   - ./figure/06_spatial_analysis/region_distances_raw.csv
#   - ./figure/06_spatial_analysis/region_distances_statistics.csv
#   - Preliminary visualization plots
#
# Author: YMS
# Date: 2024-11-17
################################################################################

# Load required packages
library(imcRtools)
library(SpatialExperiment)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)

# Set parameters
SPATIAL_PARAMS <- list(
    groups_to_compare = c("Wt", "Ko"),   # Groups for comparison
    min_samples       = 2,                # Minimum samples per group for testing
    seed              = 900223
)

# Create output directory
if (!dir.exists("./figure/06_spatial_analysis")) {
    dir.create("./figure/06_spatial_analysis", recursive = TRUE)
}

################################################################################
# 1. Load Filtered Data
################################################################################

spe <- readRDS("./data/spe_filtered.rds")

# Filter to comparison groups only
spe_region <- spe[, spe$group %in% SPATIAL_PARAMS$groups_to_compare]

cat(sprintf("Total cells: %d\n", ncol(spe_region)))
cat(sprintf("Groups: %s\n", paste(unique(spe_region$group), collapse = ", ")))

# Check if regions are defined
if (!"region" %in% colnames(colData(spe_region))) {
    stop("Column 'region' not found. Run spatial regionalization first (e.g., lisaClust).")
}

cat(sprintf("Regions: %s\n\n", 
            paste(sort(unique(na.omit(spe_region$region))), collapse = ", ")))

################################################################################
# 2. Define Cell Type Categories
################################################################################

# NOTE: Adjust these categories based on your cell type annotations
all_celltypes <- unique(spe_region$cluster_celltype)

# Epithelial cells
epithelial_types <- all_celltypes[grepl("Epithelial cell", all_celltypes)]

# Immune T cells
immune_types <- c("CTL", "CD4 T cell")
immune_types <- immune_types[immune_types %in% all_celltypes]

# Myeloid cells (if you have myeloid subclusters as SC1, SC2, etc.)
# Otherwise, use generic "Myeloid cell" annotation
if (any(grepl("^SC[0-9]+$", all_celltypes))) {
    myeloid_types <- all_celltypes[grepl("^SC[0-9]+$", all_celltypes)]
} else {
    myeloid_types <- all_celltypes[grepl("Myeloid|Monocyte|Macrophage", all_celltypes)]
}
myeloid_types <- sort(myeloid_types)

cat("Cell Type Categories:\n")
cat("=====================\n\n")
cat(sprintf("Epithelial types (%d):\n", length(epithelial_types)))
for (et in epithelial_types) cat("  *", et, "\n")
cat(sprintf("\nImmune T cell types (%d):\n", length(immune_types)))
for (it in immune_types) cat("  *", it, "\n")
cat(sprintf("\nMyeloid types (%d):\n", length(myeloid_types)))
for (mt in head(myeloid_types, 10)) cat("  *", mt, "\n")
if (length(myeloid_types) > 10) cat("  * ... and", length(myeloid_types) - 10, "more\n")

################################################################################
# 3. Define Analysis Combinations
################################################################################

# IMPORTANT: Four categories of cell-cell interactions
analysis_combinations <- list()

# Category 1: Immune → Epithelial
for (immune in immune_types) {
    for (epi in epithelial_types) {
        analysis_combinations[[length(analysis_combinations) + 1]] <- 
            list(from = immune, to = epi, category = "Immune->Epithelial")
    }
}

# Category 2: Immune → Myeloid
for (immune in immune_types) {
    for (myeloid in myeloid_types) {
        analysis_combinations[[length(analysis_combinations) + 1]] <- 
            list(from = immune, to = myeloid, category = "Immune->Myeloid")
    }
}

# Category 3: Myeloid → Immune
for (myeloid in myeloid_types) {
    for (immune in immune_types) {
        analysis_combinations[[length(analysis_combinations) + 1]] <- 
            list(from = myeloid, to = immune, category = "Myeloid->Immune")
    }
}

# Category 4: Myeloid → Epithelial
for (myeloid in myeloid_types) {
    for (epi in epithelial_types) {
        analysis_combinations[[length(analysis_combinations) + 1]] <- 
            list(from = myeloid, to = epi, category = "Myeloid->Epithelial")
    }
}

cat(sprintf("\n\nTotal analysis combinations: %d\n", length(analysis_combinations)))

################################################################################
# 4. Region-Based Distance Calculation
################################################################################

# NOTE: This may take several minutes to hours depending on data size
# Progress is tracked internally

regions <- sort(unique(na.omit(spe_region$region)))
within_results <- list()

cat("\nStarting distance calculations...\n")
start_time <- Sys.time()

for (reg in regions) {
    
    region_cells <- which(spe_region$region == reg)
    available_types <- unique(spe_region$cluster_celltype[region_cells])
    
    for (combo in analysis_combinations) {
        from_type <- combo$from
        to_type <- combo$to
        category <- combo$category
        
        # Check if both cell types exist in this region
        if (!from_type %in% available_types || !to_type %in% available_types) {
            next
        }
        
        # Get cell indices
        to_cells_idx <- which(
            spe_region$region == reg & 
            spe_region$cluster_celltype == to_type
        )
        from_cells_idx <- which(
            spe_region$region == reg & 
            spe_region$cluster_celltype == from_type
        )
        
        # Skip if too few cells
        if (length(to_cells_idx) < 5 || length(from_cells_idx) < 5) {
            next
        }
        
        tryCatch({
            # Create logical vector for target cells
            to_cells_logical <- rep(FALSE, ncol(spe_region))
            to_cells_logical[to_cells_idx] <- TRUE
            
            # Calculate minimum distances
            temp_spe <- minDistToCells(
                spe_region,
                x_cells = to_cells_logical,
                img_id = "sample_id",
                name = "temp_dist"
            )
            
            # Extract and summarize distances
            dist_data <- colData(temp_spe)[from_cells_idx, ] %>%
                as_tibble() %>%
                filter(!is.na(temp_dist)) %>%
                group_by(sample_id, group) %>%
                summarize(
                    n_cells = n(),
                    mean_dist = mean(temp_dist, na.rm = TRUE),
                    median_dist = median(temp_dist, na.rm = TRUE),
                    sd_dist = sd(temp_dist, na.rm = TRUE),
                    q25 = quantile(temp_dist, 0.25, na.rm = TRUE),
                    q75 = quantile(temp_dist, 0.75, na.rm = TRUE),
                    .groups = "drop"
                ) %>%
                mutate(
                    region = reg,
                    from_celltype = from_type,
                    to_celltype = to_type,
                    category = category
                )
            
            # Store results
            if (nrow(dist_data) > 0) {
                key <- paste(reg, from_type, to_type, sep = "___")
                within_results[[key]] <- dist_data
            }
        }, error = function(e) {
            NULL
        })
    }
}

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat(sprintf("\nDistance calculation complete: %.1f minutes\n", total_time))
cat(sprintf("Successfully calculated: %d combinations\n", length(within_results)))

################################################################################
# 5. Combine and Save Raw Results
################################################################################

if (length(within_results) > 0) {
    within_all <- bind_rows(within_results)
    
    write_csv(
        within_all,
        "./figure/06_spatial_analysis/region_distances_raw.csv"
    )
    
    cat(sprintf("Total distance measurements: %d\n", nrow(within_all)))
    
} else {
    stop("No distance results calculated. Check cell type definitions and region assignments.")
}

################################################################################
# 6. Statistical Testing
################################################################################

# Determine valid pairs (minimum sample size per group)
sample_counts <- within_all %>%
    group_by(region, from_celltype, to_celltype, group) %>%
    summarize(n_samples = n(), .groups = "drop") %>%
    pivot_wider(
        names_from = group,
        values_from = n_samples,
        values_fill = 0,
        names_prefix = "n_"
    )

# Filter for sufficient sample sizes
valid_pairs <- sample_counts %>%
    filter(if_all(starts_with("n_"), ~ . >= SPATIAL_PARAMS$min_samples))

cat(sprintf("\nValid pairs for testing: %d out of %d\n",
            nrow(valid_pairs), nrow(sample_counts)))

if (nrow(valid_pairs) > 0) {
    
    # Perform Wilcoxon test
    stats_results <- within_all %>%
        semi_join(valid_pairs, by = c("region", "from_celltype", "to_celltype")) %>%
        group_by(region, from_celltype, to_celltype, category) %>%
        wilcox_test(mean_dist ~ group) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance("p.adj") %>%
        arrange(p.adj)
    
    # Save statistics
    write_csv(
        stats_results,
        "./figure/06_spatial_analysis/region_distances_statistics.csv"
    )
    
    # Print top results
    cat("\nTop 20 Most Significant Results:\n")
    print(head(stats_results, 20), n = 20)
    
    # Summary of significant results
    sig_results <- stats_results %>% filter(p.adj < 0.05)
    
    if (nrow(sig_results) > 0) {
        cat(sprintf("\n\nSignificant pairs (p.adj < 0.05): %d\n", nrow(sig_results)))
        
        sig_by_category <- sig_results %>%
            group_by(category) %>%
            summarize(n_sig = n(), .groups = "drop") %>%
            arrange(desc(n_sig))
        
        cat("\nBy Category:\n")
        for (i in 1:nrow(sig_by_category)) {
            cat(sprintf("  * %s: %d\n",
                        sig_by_category$category[i],
                        sig_by_category$n_sig[i]))
        }
    } else {
        cat("\nNo significant results after FDR correction.\n")
        cat("Consider examining effect sizes in Step 07.\n")
    }
    
} else {
    cat("\nInsufficient sample sizes for statistical testing.\n")
}

################################################################################
# 7. Preliminary Visualization
################################################################################

# Category-wise comparison
if (nrow(within_all) > 0) {
    
    category_summary <- within_all %>%
        group_by(category, group) %>%
        summarize(
            n = n(),
            mean_dist = mean(mean_dist, na.rm = TRUE),
            se_dist = sd(mean_dist, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
        )
    
    p_category <- ggplot(
        category_summary,
        aes(x = category, y = mean_dist, fill = group)
    ) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        geom_errorbar(
            aes(ymin = mean_dist - se_dist, ymax = mean_dist + se_dist),
            position = position_dodge(0.9),
            width = 0.3
        ) +
        scale_fill_manual(
            values = metadata(spe)$color_vectors$group[SPATIAL_PARAMS$groups_to_compare]
        ) +
        theme_classic(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
        labs(
            x = "",
            y = "Mean Distance (μm)",
            title = "Cell-Cell Distance by Category",
            fill = "Group"
        )
    
    ggsave(
        filename = "./figure/06_spatial_analysis/01_category_summary.tiff",
        plot = p_category,
        width = 11,
        height = 7,
        dpi = 400,
        compression = "lzw"
    )
}

# Top pairs boxplot (if significant results exist)
if (exists("sig_results") && nrow(sig_results) > 0) {
    
    n_plots <- min(12, nrow(sig_results))
    plot_list <- list()
    
    for (i in 1:n_plots) {
        data_sub <- within_all %>%
            filter(
                region == sig_results$region[i],
                from_celltype == sig_results$from_celltype[i],
                to_celltype == sig_results$to_celltype[i]
            )
        
        p <- ggplot(data_sub, aes(x = group, y = mean_dist, fill = group)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.7) +
            geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
            scale_fill_manual(
                values = metadata(spe)$color_vectors$group[SPATIAL_PARAMS$groups_to_compare]
            ) +
            theme_classic(base_size = 9) +
            labs(
                title = sprintf("%s\n%s -> %s (p=%.1e)",
                                sig_results$region[i],
                                sig_results$from_celltype[i],
                                sig_results$to_celltype[i],
                                sig_results$p.adj[i]),
                x = "",
                y = "Distance (μm)"
            ) +
            stat_compare_means(
                method = "wilcox.test",
                label = "p.format",
                size = 2.5
            ) +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 8, face = "bold")
            )
        
        plot_list[[i]] <- p
    }
    
    p_combined <- wrap_plots(plot_list, ncol = 4)
    
    ggsave(
        filename = "./figure/06_spatial_analysis/02_significant_pairs_boxplot.tiff",
        plot = p_combined,
        width = 16,
        height = 12,
        dpi = 400,
        compression = "lzw"
    )
}

################################################################################
# 8. Clean Up
################################################################################

rm(spe_region, within_results, temp_spe)
gc()

# NOTE: Next step is 07_advanced_visualization.R
