################################################################################
# Step 02: Data Preparation and Batch Integration
#
# Purpose: Load metadata, panel info, and integrate batch-specific CyTOF data
# Dataset: NMIBC bladder cancer CyTOF
#
# Input:
#   - ./pbmc_metadata.xlsx (sample metadata with batch info)
#   - ./pbmc_panel.xlsx (panel information)
#   - ./data/sample_*.fcs (debarcoded FCS files)
#
# Output:
#   - ./data/sce_integrated.RData (integrated SCE object)
#   - ./results/02_preparation/*.png (QC plots)
#
# Author: YMS
# Date: 2025
################################################################################

library(readxl)
library(dplyr)
library(stringr)
library(flowCore)
library(CATALYST)
library(ggplot2)

# Create output directories
if (!dir.exists("./results/02_preparation")) {
    dir.create("./results/02_preparation", recursive = TRUE)
}

################################################################################
# 1. Load Metadata and Panel Information
################################################################################

# NOTE: Metadata includes sample information, clinical data, and batch identifiers
# Two batches exist:
#   - Batch A (250904): Y89 used as barcode channel
#   - Batch B (0910/0913): Y89 used as marker channel

md_all <- read_xlsx("./pbmc_metadata.xlsx") %>%
    mutate(
        recurrence = factor(
            recurrence,
            levels = c("control", "recurrence", "non-recurrence")
        ),
        condition = factor(
            condition,
            levels = c("ctrl", "UR", "RE", "RC", "DX")
        ),
        batch = case_when(
            str_starts(file_name, "250904") ~ "A_Y89_barcode",
            str_starts(file_name, "0910")   ~ "B_Y89_marker",
            str_starts(file_name, "0913")   ~ "B_Y89_marker",
            TRUE ~ "B_Y89_marker"
        ),
        condition = factor(
            ifelse(
                is.na(condition) & str_starts(sample_id, "Control"),
                "ctrl",
                as.character(condition)
            ),
            levels = c("ctrl", "UR", "RE", "RC", "DX")
        )
    )

# Split metadata by batch
md_A <- dplyr::filter(md_all, batch == "A_Y89_barcode")
md_B <- dplyr::filter(md_all, batch == "B_Y89_marker")

# Load panel information
panel_base <- read_xlsx("./pbmc_panel.xlsx") %>%
    mutate(
        marker_class = tolower(marker_class),
        use_channel  = as.logical(use_channel)
    )

################################################################################
# 2. Prepare Batch-Specific Panels
################################################################################

# NOTE: Barcode channels differ between batches
# Must handle Y89 appropriately for each batch

bc_chs_all <- c("Y89Di", "Pt194Di", "Pt195Di", "Pt196Di", "Pt198Di", "Cd113Di")

# Function to mark barcode channels
make_panel_bc <- function(panel, bc_chs) {
    panel %>%
        mutate(
            antigen = ifelse(
                fcs_colname %in% bc_chs,
                paste0("BC_CD45_", fcs_colname),
                antigen
            ),
            marker_class = ifelse(
                fcs_colname %in% bc_chs,
                "none",
                marker_class
            ),
            use_channel = ifelse(
                fcs_colname %in% bc_chs,
                FALSE,
                use_channel
            )
        )
}

# Batch A: All 6 channels as barcodes (including Y89)
panel_A <- make_panel_bc(panel_base, bc_chs_all)

# Batch B: 5 channels as barcodes (excluding Y89, kept as marker)
panel_B <- make_panel_bc(panel_base, setdiff(bc_chs_all, "Y89Di"))

################################################################################
# 3. Load FCS Files and Create SCE Objects
################################################################################

# Function to create flowSet
make_flowset <- function(md, panel, fcs_path = "./data") {
    chs  <- panel$fcs_colname[!is.na(panel$fcs_colname)]
    patt <- paste(chs, collapse = "|")
    
    read.flowSet(
        path = fcs_path,
        transformation = FALSE,
        truncate_max_range = FALSE,
        files = md$file_name,
        column.pattern = patt
    )
}

# Function to create SCE with robust feature matching
make_sce <- function(fs, panel, md) {
    # NOTE: Match only channels that exist in both panel and FCS files
    existing <- intersect(panel$fcs_colname, colnames(fs[[1]]))
    panel2   <- dplyr::filter(panel, fcs_colname %in% existing)
    feats    <- panel2$fcs_colname
    
    prepData(
        x          = fs,
        panel      = panel2,
        md         = md,
        features   = feats,
        md_cols    = list(
            file = "file_name",
            id = "sample_id",
            factors = intersect(
                c("recurrence", "condition", "group", "status"),
                names(md)
            )
        ),
        panel_cols = list(
            channel = "fcs_colname",
            antigen = "antigen",
            marker_class = "marker_class"
        )
    )
}

# Create batch-specific SCE objects
fs_A  <- make_flowset(md_A, panel_A)
sce_A <- make_sce(fs_A, panel_A, md_A)

fs_B  <- make_flowset(md_B, panel_B)
sce_B <- make_sce(fs_B, panel_B, md_B)

################################################################################
# 4. Merge Batches on Common Markers
################################################################################

# NOTE: Only markers with use_channel == TRUE are used for analysis
feat_A <- rownames(sce_A)[rowData(sce_A)$use_channel]
feat_B <- rownames(sce_B)[rowData(sce_B)$use_channel]

# Find common markers between batches
shared <- intersect(feat_A, feat_B)

# Merge batches keeping only shared markers
sce <- cbind(sce_A[shared, ], sce_B[shared, ])

# Print integration summary
cat("=== Batch Integration Summary ===\n")
cat("Batch A markers:", length(feat_A), "\n")
cat("Batch B markers:", length(feat_B), "\n")
cat("Shared markers:", length(shared), "\n")
cat("Final SCE markers:", nrow(sce), "\n")
cat("Final SCE cells:", ncol(sce), "\n\n")

# Print cell counts per sample
print(n_cells(sce))

################################################################################
# 5. Add Batch Information to SCE
################################################################################

# Add detailed batch and Y89 role information
md_all <- md_all %>%
    mutate(
        batch_run = sub("_sample.*", "", file_name),
        y89_role  = ifelse(
            grepl("^250904", file_name),
            "Y89_barcode",
            "Y89_marker"
        )
    )

# Match and add to SCE colData
key <- md_all[, c("sample_id", "batch_run", "y89_role")]
idx <- match(colData(sce)$sample_id, key$sample_id)

colData(sce)$batch <- factor(
    key$batch_run[idx],
    levels = sort(unique(key$batch_run))
)
colData(sce)$y89_role <- factor(
    key$y89_role[idx],
    levels = c("Y89_barcode", "Y89_marker")
)

# Check batch distribution
cat("=== Batch Distribution ===\n")
print(table(colData(sce)$batch))
cat("\n")

cat("=== Batch vs Condition Distribution ===\n")
print(table(colData(sce)$batch, colData(sce)$condition))

################################################################################
# 6. Quality Control Plots
################################################################################

# Plot number of cells per sample
nrs_plot <- plotNRS(sce, features = "type", color_by = "group")
ggsave(
    filename = "./results/02_preparation/cells_per_sample.png",
    plot = nrs_plot,
    width = 10,
    height = 6,
    dpi = 300
)

################################################################################
# 7. Filter Samples (Optional)
################################################################################

# NOTE: Exclude low-quality or problematic samples if needed
# Example filtering (adjust sample IDs as needed):

excluded_samples <- c(
    "NMIBC_04", "NMIBC_14", "NMIBC_17", "NMIBC_18",
    "NMIBC_19", "NMIBC_24", "NMIBC_25", "NMIBC_26",
    "NMIBC_27", "NMIBC_42"
)

# Save original SCE
sce_original <- sce

# Apply filtering if needed
if (length(excluded_samples) > 0) {
    sce <- filterSCE(sce, !(sample_id %in% excluded_samples))
    cat("\n=== Sample Filtering ===\n")
    cat("Original samples:", ncol(sce_original), "\n")
    cat("Filtered samples:", ncol(sce), "\n")
    cat("Excluded:", length(excluded_samples), "samples\n")
}

################################################################################
# 8. Save Integrated Data
################################################################################

save(sce, file = "./data/sce_integrated.RData")
save(sce_original, file = "./data/sce_original.RData")

# Clean up
rm(fs_A, fs_B, sce_A, sce_B)
gc()

################################################################################
# End of Step 02
################################################################################
