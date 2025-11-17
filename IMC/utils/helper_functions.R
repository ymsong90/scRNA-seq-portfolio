################################################################################
# IMC Analysis Pipeline
# Helper Functions
# 
# Description:
#   - Utility functions used across multiple pipeline steps
#   - Color palette generators
#   - Data transformation helpers
#   - Plotting utilities
#   - Quality control functions
#
# Usage:
#   source("IMC/utils/helper_functions.R")
#
# Author: YMS
# Date: 2024-11-17
################################################################################

################################################################################
# 1. Color Palette Functions
################################################################################

#' Generate color palette for clusters
#' 
#' @param n Number of colors needed
#' @return Character vector of hex colors
#' 
make_cluster_colors <- function(n) {
    pals <- c("Set3", "Paired", "Set1", "Set2", "Accent", "Dark2", "Pastel1", "Pastel2")
    cols <- unlist(lapply(pals, function(p) {
        k <- RColorBrewer::brewer.pal.info[p, "maxcolors"]
        RColorBrewer::brewer.pal(as.integer(k), p)
    }))
    cols <- unique(cols)
    if (length(cols) < n) {
        cols <- c(cols, scales::hue_pal()(n - length(cols)))
    }
    cols[seq_len(n)]
}

#' Get group colors (standard palette)
#' 
#' @param groups Character vector of group names
#' @return Named vector of colors
#' 
get_group_colors <- function(groups) {
    base_colors <- c(
        Wt = "#4e79a7",
        Ko = "#e15759",
        Panc = "#59a14f",
        WT = "#4e79a7",
        KO = "#e15759"
    )
    
    # Add colors for groups not in base palette
    missing_groups <- setdiff(groups, names(base_colors))
    if (length(missing_groups) > 0) {
        extra_colors <- setNames(
            scales::hue_pal()(length(missing_groups)),
            missing_groups
        )
        base_colors <- c(base_colors, extra_colors)
    }
    
    base_colors[groups]
}

################################################################################
# 2. Data Transformation Functions
################################################################################

#' Asinh transformation for IMC data
#' 
#' @param x Numeric vector or matrix of counts
#' @param cofactor Transformation cofactor (default: 5)
#' @return Transformed values
#' 
asinh_transform <- function(x, cofactor = 5) {
    asinh(x / cofactor)
}

#' Inverse asinh transformation
#' 
#' @param x Transformed values
#' @param cofactor Transformation cofactor (default: 5)
#' @return Original scale values
#' 
inv_asinh_transform <- function(x, cofactor = 5) {
    sinh(x) * cofactor
}

#' Z-score normalization with clipping
#' 
#' @param mat Matrix (features x samples)
#' @param clip Clipping threshold (default: 3)
#' @param byrow Normalize by row (default: TRUE)
#' @return Z-score normalized matrix
#' 
zscore_normalize <- function(mat, clip = 3, byrow = TRUE) {
    if (byrow) {
        mat_z <- t(scale(t(mat)))
    } else {
        mat_z <- scale(mat)
    }
    mat_z[is.na(mat_z)] <- 0
    mat_z <- pmin(pmax(mat_z, -clip), clip)
    mat_z
}

################################################################################
# 3. Quality Control Functions
################################################################################

#' Calculate percentage of cells expressing a marker
#' 
#' @param spe SpatialExperiment object
#' @param marker Marker name
#' @param threshold Expression threshold (default: 0.5 on exprs scale)
#' @return Data frame with percentages per sample
#' 
calculate_marker_percentage <- function(spe, marker, threshold = 0.5) {
    if (!marker %in% rownames(spe)) {
        stop(sprintf("Marker '%s' not found in SPE object", marker))
    }
    
    expr <- assay(spe, "exprs")[marker, ]
    positive <- expr > threshold
    
    result <- colData(spe) %>%
        as.data.frame() %>%
        mutate(positive = positive) %>%
        group_by(sample_id, group) %>%
        summarize(
            n_total = n(),
            n_positive = sum(positive),
            pct_positive = 100 * n_positive / n_total,
            .groups = "drop"
        )
    
    result
}

#' Calculate neighborhood enrichment score
#' 
#' @param spe SpatialExperiment object
#' @param from_celltype Source cell type
#' @param to_celltype Target cell type
#' @return Enrichment score (observed/expected ratio)
#' 
calculate_neighborhood_enrichment <- function(spe, from_celltype, to_celltype) {
    if (!"neighborhood" %in% names(colPairs(spe))) {
        stop("Neighborhood information not found. Run buildSpatialGraph first.")
    }
    
    pairs <- colPair(spe, "neighborhood")
    from_idx <- which(spe$cluster_celltype == from_celltype)
    to_idx <- which(spe$cluster_celltype == to_celltype)
    
    # Observed interactions
    observed <- sum(
        from(pairs) %in% from_idx &
        to(pairs) %in% to_idx
    )
    
    # Expected interactions (random)
    total_pairs <- length(pairs)
    p_from <- length(from_idx) / ncol(spe)
    p_to <- length(to_idx) / ncol(spe)
    expected <- total_pairs * p_from * p_to
    
    # Enrichment score
    enrichment <- observed / expected
    
    list(
        observed = observed,
        expected = expected,
        enrichment = enrichment
    )
}

################################################################################
# 4. Plotting Utilities
################################################################################

#' Create publication-quality theme
#' 
#' @param base_size Base font size (default: 12)
#' @return ggplot2 theme object
#' 
theme_publication <- function(base_size = 12) {
    ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
            axis.line = ggplot2::element_line(size = 0.5, color = "black"),
            axis.text = ggplot2::element_text(color = "black", size = base_size * 0.9),
            axis.title = ggplot2::element_text(face = "bold", size = base_size * 1.1),
            legend.title = ggplot2::element_text(face = "bold", size = base_size),
            legend.text = ggplot2::element_text(size = base_size * 0.9),
            strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
            strip.text = ggplot2::element_text(face = "bold", size = base_size),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.2),
            plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30")
        )
}

#' Save plot with consistent settings
#' 
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Width in inches (default: 10)
#' @param height Height in inches (default: 8)
#' @param dpi Resolution (default: 400)
#' 
save_publication_plot <- function(plot, filename, width = 10, height = 8, dpi = 400) {
    # Ensure directory exists
    dir_name <- dirname(filename)
    if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Save as TIFF with LZW compression
    tiff_file <- sub("\\.pdf$|\\.png$|\\.tiff$", ".tiff", filename)
    ggplot2::ggsave(
        filename = tiff_file,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        compression = "lzw"
    )
    
    # Also save as PDF (vector format)
    pdf_file <- sub("\\.tiff$|\\.png$|\\.pdf$", ".pdf", filename)
    ggplot2::ggsave(
        filename = pdf_file,
        plot = plot,
        width = width,
        height = height
    )
    
    invisible(list(tiff = tiff_file, pdf = pdf_file))
}

################################################################################
# 5. Data Export Functions
################################################################################

#' Export cell type proportions
#' 
#' @param spe SpatialExperiment object
#' @param group_by Grouping variable (default: "sample_id")
#' @param filename Output CSV filename
#' 
export_proportions <- function(spe, group_by = "sample_id", filename = NULL) {
    prop_data <- colData(spe) %>%
        as.data.frame() %>%
        group_by(across(all_of(c(group_by, "cluster_celltype")))) %>%
        summarize(n = n(), .groups = "drop") %>%
        group_by(across(all_of(group_by)))
