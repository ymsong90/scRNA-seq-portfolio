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
#' IMPORTANT: RColorBrewer의 모든 팔레트를 결합하여 충분한 색상 생성
#' 클러스터 수가 많을 때 유용 (35개 이상)
#' 
#' @param n Number of colors needed
#' @return Character vector of hex colors
#' 
make_cluster_colors <- function(n) {
    # NOTE: 사용 가능한 모든 ColorBrewer 팔레트 활용
    pals <- c("Set3", "Paired", "Set1", "Set2", "Accent", "Dark2", "Pastel1", "Pastel2")
    
    cols <- unlist(lapply(pals, function(p) {
        k <- RColorBrewer::brewer.pal.info[p, "maxcolors"]
        RColorBrewer::brewer.pal(as.integer(k), p)
    }))
    
    # 중복 제거
    cols <- unique(cols)
    
    # IMPORTANT: 필요한 색상이 더 많으면 hue_pal()로 추가 생성
    if (length(cols) < n) {
        cols <- c(cols, scales::hue_pal()(n - length(cols)))
    }
    
    cols[seq_len(n)]
}

#' Get group colors (standard palette)
#' 
#' IMPORTANT: 일관된 그룹 색상 사용 (모든 플롯에서 동일)
#' Wt = 파랑, Ko = 빨강, Panc = 녹색
#' 
#' @param groups Character vector of group names
#' @return Named vector of colors
#' 
get_group_colors <- function(groups) {
    # 기본 그룹 색상 (대소문자 모두 지원)
    base_colors <- c(
        Wt = "#4e79a7",
        Ko = "#e15759",
        Panc = "#59a14f",
        WT = "#4e79a7",
        KO = "#e15759"
    )
    
    # NOTE: 기본 팔레트에 없는 그룹은 자동으로 색상 할당
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
#' IMPORTANT: IMC 데이터의 표준 변환 방법
#' CyTOF/IMC에서 cofactor=5가 일반적으로 사용됨
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
#' NOTE: 변환된 값을 원래 스케일로 복원할 때 사용
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
#' IMPORTANT: 히트맵 생성 시 극단값 제거를 위한 clipping 적용
#' 기본 clip=3은 ±3 표준편차 이상 값을 제한
#' 
#' @param mat Matrix (features x samples)
#' @param clip Clipping threshold (default: 3)
#' @param byrow Normalize by row (default: TRUE)
#' @return Z-score normalized matrix
#' 
zscore_normalize <- function(mat, clip = 3, byrow = TRUE) {
    # Row-wise 또는 column-wise normalization
    if (byrow) {
        mat_z <- t(scale(t(mat)))
    } else {
        mat_z <- scale(mat)
    }
    
    # NA를 0으로 변환 (표준편차가 0인 경우)
    mat_z[is.na(mat_z)] <- 0
    
    # IMPORTANT: Clip extreme values for better visualization
    mat_z <- pmin(pmax(mat_z, -clip), clip)
    mat_z
}

################################################################################
# 3. Quality Control Functions
################################################################################

#' Calculate percentage of cells expressing a marker
#' 
#' NOTE: 특정 마커를 발현하는 세포의 비율 계산
#' QC 및 마커 검증에 유용
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
    
    # IMPORTANT: exprs assay 사용 (asinh 변환된 값)
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
#' IMPORTANT: 두 cell type이 이웃으로 있을 확률의 enrichment 계산
#' >1 = enriched (예상보다 많이 이웃), <1 = depleted
#' 
#' @param spe SpatialExperiment object
#' @param from_celltype Source cell type
#' @param to_celltype Target cell type
#' @return Enrichment score (observed/expected ratio)
#' 
calculate_neighborhood_enrichment <- function(spe, from_celltype, to_celltype) {
    # NOTE: buildSpatialGraph() 함수로 사전에 neighborhood 생성 필요
    if (!"neighborhood" %in% names(colPairs(spe))) {
        stop("Neighborhood information not found. Run buildSpatialGraph first.")
    }
    
    pairs <- colPair(spe, "neighborhood")
    from_idx <- which(spe$cluster_celltype == from_celltype)
    to_idx <- which(spe$cluster_celltype == to_celltype)
    
    # Observed interactions: 실제 이웃 관계 수
    observed <- sum(
        from(pairs) %in% from_idx &
        to(pairs) %in% to_idx
    )
    
    # Expected interactions: 랜덤 분포일 때 예상되는 수
    total_pairs <- length(pairs)
    p_from <- length(from_idx) / ncol(spe)
    p_to <- length(to_idx) / ncol(spe)
    expected <- total_pairs * p_from * p_to
    
    # IMPORTANT: Enrichment score = observed / expected
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
#' IMPORTANT: 저널 출판을 위한 표준 테마
#' 깔끔한 배경, 명확한 축, 적절한 폰트 크기
#' 
#' @param base_size Base font size (default: 12)
#' @return ggplot2 theme object
#' 
theme_publication <- function(base_size = 12) {
    ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
            # Axis
            axis.line = ggplot2::element_line(size = 0.5, color = "black"),
            axis.text = ggplot2::element_text(color = "black", size = base_size * 0.9),
            axis.title = ggplot2::element_text(face = "bold", size = base_size * 1.1),
            
            # Legend
            legend.title = ggplot2::element_text(face = "bold", size = base_size),
            legend.text = ggplot2::element_text(size = base_size * 0.9),
            
            # Facet
            strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
            strip.text = ggplot2::element_text(face = "bold", size = base_size),
            
            # Grid (removed for cleaner look)
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            
            # Title
            plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.2),
            plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30")
        )
}

#' Save plot with consistent settings
#' 
#' IMPORTANT: TIFF (고해상도 래스터) + PDF (벡터) 동시 저장
#' 대부분의 저널이 요구하는 형식
#' 
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Width in inches (default: 10)
#' @param height Height in inches (default: 8)
#' @param dpi Resolution (default: 400)
#' @return List of saved file paths
#' 
save_publication_plot <- function(plot, filename, width = 10, height = 8, dpi = 400) {
    # Ensure directory exists
    dir_name <- dirname(filename)
    if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)
    }
    
    # IMPORTANT: TIFF with LZW compression (고해상도, 무손실 압축)
    tiff_file <- sub("\\.pdf$|\\.png$|\\.tiff$", ".tiff", filename)
    ggplot2::ggsave(
        filename = tiff_file,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        compression = "lzw"
    )
    
    # PDF for vector graphics (확대해도 깨끗)
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
#' NOTE: 샘플별 또는 그룹별 세포 타입 비율 계산 및 저장
#' 
#' @param spe SpatialExperiment object
#' @param group_by Grouping variable (default: "sample_id")
#' @param filename Output CSV filename (optional)
#' @return Data frame with proportions
#' 
export_proportions <- function(spe, group_by = "sample_id", filename = NULL) {
    prop_data <- colData(spe) %>%
        as.data.frame() %>%
        group_by(across(all_of(c(group_by, "cluster_celltype")))) %>%
        summarize(n = n(), .groups = "drop") %>%
        group_by(across(all_of(group_by))) %>%
        mutate(
            total = sum(n),
            proportion = n / total,
            percentage = 100 * proportion
        ) %>%
        ungroup()
    
    # Save to CSV if filename provided
    if (!is.null(filename)) {
        write.csv(prop_data, filename, row.names = FALSE)
    }
    
    prop_data
}

#' Export marker expression summary
#' 
#' NOTE: 세포 타입별 마커 발현 통계 (평균, 중앙값, 표준편차)
#' 
#' @param spe SpatialExperiment object
#' @param markers Character vector of marker names (NULL = all use_channel markers)
#' @param group_by Grouping variable (default: "cluster_celltype")
#' @param filename Output CSV filename (optional)
#' @return Data frame with expression summaries
#' 
export_marker_summary <- function(spe, markers = NULL, 
                                 group_by = "cluster_celltype",
                                 filename = NULL) {
    # IMPORTANT: 기본적으로 분석에 사용된 마커만 사용
    if (is.null(markers)) {
        markers <- rownames(spe)[rowData(spe)$use_channel]
    }
    
    expr_mat <- assay(spe, "exprs")[markers, , drop = FALSE]
    
    summary_data <- colData(spe) %>%
        as.data.frame() %>%
        select(all_of(group_by)) %>%
        bind_cols(as.data.frame(t(expr_mat))) %>%
        pivot_longer(
            cols = -all_of(group_by),
            names_to = "marker",
            values_to = "expression"
        ) %>%
        group_by(across(all_of(c(group_by, "marker")))) %>%
        summarize(
            n = n(),
            mean_expr = mean(expression, na.rm = TRUE),
            median_expr = median(expression, na.rm = TRUE),
            sd_expr = sd(expression, na.rm = TRUE),
            .groups = "drop"
        )
    
    if (!is.null(filename)) {
        write.csv(summary_data, filename, row.names = FALSE)
    }
    
    summary_data
}

################################################################################
# 6. Spatial Analysis Helpers
################################################################################

#' Calculate pairwise distances between cell types
#' 
#' IMPORTANT: 각 "from" 세포에서 가장 가까운 "to" 세포까지의 거리 계산
#' Euclidean distance 사용 (공간 좌표 기반)
#' 
#' @param spe SpatialExperiment object
#' @param from_celltype Source cell type
#' @param to_celltype Target cell type
#' @return Vector of minimum distances for each "from" cell
#' 
calculate_pairwise_distances <- function(spe, from_celltype, to_celltype) {
    from_idx <- which(spe$cluster_celltype == from_celltype)
    to_idx <- which(spe$cluster_celltype == to_celltype)
    
    # NOTE: 세포가 없으면 빈 벡터 반환
    if (length(from_idx) == 0 || length(to_idx) == 0) {
        return(numeric(0))
    }
    
    # Extract spatial coordinates
    from_coords <- spatialCoords(spe)[from_idx, , drop = FALSE]
    to_coords <- spatialCoords(spe)[to_idx, , drop = FALSE]
    
    # IMPORTANT: 전체 거리 행렬 계산 (computationally expensive for large datasets)
    dist_mat <- as.matrix(dist(rbind(from_coords, to_coords)))
    n_from <- nrow(from_coords)
    
    # Extract distances from "from" cells to "to" cells
    distances <- dist_mat[1:n_from, (n_from + 1):ncol(dist_mat), drop = FALSE]
    
    # Return minimum distance for each "from" cell
    apply(distances, 1, min)
}

################################################################################
# 7. Validation Functions
################################################################################

#' Check if SPE object has required components
#' 
#' NOTE: 파이프라인 진행 전 SPE 객체 유효성 검증
#' 
#' @param spe SpatialExperiment object
#' @param required Vector of required assay names
#' @return Logical indicating if all required components exist
#' 
validate_spe <- function(spe, required = c("exprs", "counts")) {
    checks <- list(
        assays = all(required %in% assayNames(spe)),
        coords = !is.null(spatialCoords(spe)),
        coldata = ncol(colData(spe)) > 0,
        rowdata = ncol(rowData(spe)) > 0
    )
    
    all(unlist(checks))
}

#' Check for batch effects
#' 
#' IMPORTANT: Silhouette score로 배치 효과 정도 측정
#' Score가 높을수록 배치 간 분리가 명확 (배치 효과가 큼)
#' 
#' @param spe SpatialExperiment object
#' @param batch_var Batch variable name (e.g., "sample_id")
#' @param reduction Dimensionality reduction to use (default: "PCA")
#' @return Silhouette score (higher = more separation)
#' 
check_batch_effect <- function(spe, batch_var, reduction = "PCA") {
    if (!reduction %in% names(reducedDims(spe))) {
        stop(sprintf("Reduction '%s' not found", reduction))
    }
    
    # NOTE: 처음 2개 PC만 사용 (시각화 목적)
    coords <- reducedDim(spe, reduction)[, 1:2]
    batch <- colData(spe)[[batch_var]]
    
    if (is.null(batch)) {
        stop(sprintf("Batch variable '%s' not found", batch_var))
    }
    
    # IMPORTANT: Silhouette coefficient 계산
    # -1 to 1 범위: 1에 가까울수록 잘 분리됨 (배치 효과가 큼)
    dist_mat <- dist(coords)
    sil <- cluster::silhouette(as.integer(factor(batch)), dist_mat)
    
    mean(sil[, "sil_width"])
}

################################################################################
# 8. Statistical Helper Functions
################################################################################

#' Calculate Cohen's d effect size
#' 
#' IMPORTANT: 두 그룹 간 차이의 크기를 표준화된 단위로 표현
#' 해석: 0.2=small, 0.5=medium, 0.8=large effect
#' 
#' @param x1 Numeric vector for group 1
#' @param x2 Numeric vector for group 2
#' @return Cohen's d value
#' 
cohens_d <- function(x1, x2) {
    n1 <- length(x1)
    n2 <- length(x2)
    
    # Pooled standard deviation
    pooled_sd <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))
    
    # Cohen's d
    d <- (mean(x1) - mean(x2)) / pooled_sd
    
    return(d)
}

#' Perform multiple testing correction
#' 
#' NOTE: FDR (False Discovery Rate) 조정으로 다중 검정 문제 해결
#' 
#' @param p_values Vector of p-values
#' @param method Correction method (default: "BH" for Benjamini-Hochberg)
#' @return Adjusted p-values
#' 
adjust_pvalues <- function(p_values, method = "BH") {
    p.adjust(p_values, method = method)
}

################################################################################
# 9. Spatial Regionalization Helpers
################################################################################

#' Calculate spatial autocorrelation (Moran's I)
#' 
#' IMPORTANT: 공간적 패턴 존재 여부 확인
#' Moran's I > 0: 유사한 값이 근처에 모임 (positive autocorrelation)
#' Moran's I < 0: 다른 값이 근처에 모임 (negative autocorrelation)
#' 
#' @param spe SpatialExperiment object
#' @param feature Feature name (marker or cell type)
#' @param k Number of neighbors (default: 10)
#' @return Moran's I statistic
#' 
calculate_morans_i <- function(spe, feature, k = 10) {
    # Get spatial coordinates
    coords <- spatialCoords(spe)
    
    # Get feature values
    if (feature %in% rownames(spe)) {
        values <- assay(spe, "exprs")[feature, ]
    } else if (feature %in% colnames(colData(spe))) {
        values <- colData(spe)[[feature]]
        if (!is.numeric(values)) {
            values <- as.numeric(factor(values))
        }
    } else {
        stop("Feature not found in SPE object")
    }
    
    # IMPORTANT: k-nearest neighbors graph 생성
    knn <- RANN::nn2(coords, k = k + 1)$nn.idx[, -1]
    
    # Create weights matrix
    n <- length(values)
    W <- matrix(0, n, n)
    for (i in 1:n) {
        W[i, knn[i, ]] <- 1
    }
    
    # Calculate Moran's I
    W_sum <- sum(W)
    z <- values - mean(values)
    numerator <- sum(W * outer(z, z))
    denominator <- sum(z^2)
    
    I <- (n / W_sum) * (numerator / denominator)
    
    return(I)
}

################################################################################
# End of helper_functions.R
################################################################################
