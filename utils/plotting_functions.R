################################################################################
# Utility Functions: Plotting
#
# Purpose: Common plotting functions for scRNA-seq analysis
# Author: YMS
# Date: 2025
################################################################################

library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

################################################################################
# 1. Enhanced Feature Plot
################################################################################

#' Generate high-quality feature plots
#'
#' @param seu Seurat object
#' @param features Vector of gene names
#' @param output_dir Output directory
#' @param pt_size Point size (default: 0.5)
#' @param use_custom_colors Use custom color palette
#'
plot_features_high_quality <- function(seu, features, output_dir,
                                      pt_size = 0.5,
                                      use_custom_colors = TRUE) {
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Define color palette
    if (use_custom_colors) {
        palette_full <- brewer.pal(9, "YlOrRd")
        palette_adj  <- palette_full[3:9]
    } else {
        palette_adj <- c("grey", "blue")
    }
    
    for (gene in features) {
        
        if (!gene %in% rownames(seu)) {
            warning(paste("Gene", gene, "not found. Skipping."))
            next
        }
        
        p <- FeaturePlot(
            object     = seu,
            features   = gene,
            cols       = palette_adj,
            reduction  = "umap",
            pt.size    = pt_size,
            order      = TRUE
        ) +
            labs(title = gene) +
            theme(
                plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank()
            )
        
        ggsave(
            filename = paste0(output_dir, "/", gene, "_FeaturePlot.tiff"),
            plot = p,
            device = "tiff",
            width = 10,
            height = 10,
            dpi = 600,
            compression = "lzw"
        )
    }
}

################################################################################
# 2. Custom Volcano Plot
################################################################################

#' Generate volcano plot for DEG results
#'
#' @param deg_results Data frame from FindMarkers
#' @param title Plot title
#' @param lfc_cutoff Log fold-change cutoff
#' @param padj_cutoff Adjusted p-value cutoff
#' @param label_n Number of top genes to label
#'
plot_volcano_custom <- function(deg_results, title = "Volcano Plot",
                               lfc_cutoff = 0.25, padj_cutoff = 0.05,
                               label_n = 20) {
    
    # Prepare data
    df <- deg_results %>%
        mutate(
            log10P = -log10(p_val_adj + 1e-300),
            sig_flag = case_when(
                p_val_adj < padj_cutoff & avg_log2FC > lfc_cutoff ~ "Up",
                p_val_adj < padj_cutoff & avg_log2FC < -lfc_cutoff ~ "Down",
                TRUE ~ "NS"
            ),
            sig_flag = factor(sig_flag, levels = c("Up", "Down", "NS"))
        )
    
    # Select top genes to label
    top_genes <- df %>%
        filter(sig_flag != "NS") %>%
        arrange(p_val_adj) %>%
        group_by(sig_flag) %>%
        slice_head(n = label_n) %>%
        ungroup()
    
    # Create plot
    p <- ggplot(df, aes(x = avg_log2FC, y = log10P, color = sig_flag)) +
        geom_point(size = 2, alpha = 0.6) +
        scale_color_manual(
            breaks = c("Up", "Down", "NS"),
            values = c("Up" = "#f06292", "Down" = "#311b92", "NS" = "grey50"),
            labels = c(
                paste0("Up (", sum(df$sig_flag == "Up"), ")"),
                paste0("Down (", sum(df$sig_flag == "Down"), ")"),
                "NS"
            )
        ) +
        geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                  linetype = "dashed", color = "black", alpha = 0.5) +
        geom_hline(yintercept = -log10(padj_cutoff),
                  linetype = "dashed", color = "black", alpha = 0.5) +
        geom_text_repel(
            data = top_genes,
            aes(label = gene),
            size = 3,
            max.overlaps = 50,
            show.legend = FALSE
        ) +
        labs(
            title = title,
            x = "log2 Fold Change (KO vs WT)",
            y = "-log10(Adjusted P-value)",
            color = "Regulation"
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "right",
            panel.grid.minor = element_blank()
        )
    
    return(p)
}

################################################################################
# 3. Enhanced DotPlot
################################################################################

#' Create enhanced dot plot with custom styling
#'
#' @param seu Seurat object
#' @param features Vector of genes
#' @param group_by Metadata column to group by
#' @param colors Color scheme
#'
dotplot_enhanced <- function(seu, features, group_by = NULL,
                            colors = "RdYlBu", dot_scale = 10) {
    
    if (!is.null(group_by)) {
        Idents(seu) <- group_by
    }
    
    p <- DotPlot(
        object = seu,
        features = features,
        cols = colors,
        dot.scale = dot_scale
    ) +
        RotatedAxis() +
        theme_bw(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                      colour = "black", size = 10),
            axis.text.y = element_text(colour = "black", size = 10),
            axis.ticks = element_blank(),
            legend.title = element_text(size = 12, face = "bold"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
        ) +
        labs(x = NULL, y = NULL)
    
    return(p)
}

################################################################################
# 4. Save Plot in Multiple Formats
################################################################################

#' Save plot in multiple formats
#'
#' @param plot ggplot object
#' @param filename Base filename (without extension)
#' @param width Plot width
#' @param height Plot height
#' @param formats Vector of formats
#'
save_plot_multi_format <- function(plot, filename, width = 10, height = 8,
                                   formats = c("png", "tiff"), dpi = 600) {
    
    for (fmt in formats) {
        out_file <- paste0(filename, ".", fmt)
        
        if (fmt == "tiff") {
            ggsave(
                filename = out_file,
                plot = plot,
                device = "tiff",
                width = width,
                height = height,
                dpi = dpi,
                compression = "lzw"
            )
        } else {
            ggsave(
                filename = out_file,
                plot = plot,
                device = fmt,
                width = width,
                height = height,
                dpi = dpi
            )
        }
    }
}

################################################################################
# 5. Create Cell Type Color Palette
################################################################################

#' Generate color palette for cell types
#'
#' @param cell_types Vector of cell type names
#' @param palette_name RColorBrewer palette name
#'
create_celltype_palette <- function(cell_types, palette_name = NULL) {
    
    n_types <- length(cell_types)
    
    if (is.null(palette_name)) {
        if (n_types <= 12) {
            colors <- brewer.pal(max(3, n_types), "Set3")[1:n_types]
        } else {
            pal1 <- brewer.pal(12, "Set3")
            pal2 <- brewer.pal(8, "Pastel2")
            colors <- c(pal1, pal2)[1:n_types]
        }
    } else {
        max_colors <- brewer.pal.info[palette_name, "maxcolors"]
        colors <- colorRampPalette(brewer.pal(max_colors, palette_name))(n_types)
    }
    
    names(colors) <- cell_types
    return(colors)
}

################################################################################
# End of Plotting Functions
################################################################################
