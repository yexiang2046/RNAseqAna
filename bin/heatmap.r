#!/usr/bin/env Rscript

# Heatmap of significant DEGs from edger.r output
# Input: DEG CSV + normalized_counts_CPM.csv
# Thresholds match functional_analysis.r defaults: FDR < 0.05 AND |logFC| > 1

suppressPackageStartupMessages({
  library(optparse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-d", "--deg"),
              type = "character",
              help = "Comma-separated DEG CSV files from edger.r. Significant genes are unioned across all files."),
  make_option(c("-n", "--norm_counts"),
              type = "character",
              help = "Comma-separated normalized counts CPM CSV files from edger.r (one or more normalized_counts_CPM.csv)"),
  make_option(c("-m", "--metadata"),
              type = "character", default = NULL,
              help = "Metadata file (tab-delimited, SampleId + group columns) for column annotation"),
  make_option(c("-o", "--output"),
              type = "character", default = ".",
              help = "Output directory [default: current directory]"),
  make_option(c("-f", "--fdr"),
              type = "double", default = 0.05,
              help = "FDR threshold [default: 0.05]"),
  make_option(c("-l", "--logfc"),
              type = "double", default = 1.0,
              help = "Minimum |logFC| threshold [default: 1.0]"),
  make_option("--use_logfc",
              type = "character", default = "TRUE",
              help = "Apply |logFC| threshold in addition to FDR (TRUE/FALSE) [default: TRUE]"),
  make_option("--max_genes",
              type = "integer", default = 100,
              help = "Maximum number of genes to plot, ranked by |logFC| [default: 100]"),
  make_option("--samples",
              type = "character", default = NULL,
              help = "Comma-separated sample IDs to include (default: all)"),
  make_option("--gap_samples",
              type = "character", default = NULL,
              help = paste("Comma-separated sample IDs after which to draw a column gap.",
                           "Disables column clustering. Use with --samples to fix column order.",
                           "Example: --gap_samples 'WT_rep3,KO_rep3'")),
  make_option("--prefix",
              type = "character", default = NULL,
              help = "Output file prefix (default: derived from --deg filename)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$deg))         stop("--deg is required")
if (is.null(opt$norm_counts)) stop("--norm_counts is required")

use_logfc <- !tolower(opt$use_logfc) %in% c("false", "0", "no", "n")

dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

deg_files <- trimws(strsplit(opt$deg, ",")[[1]])

prefix <- opt$prefix
if (is.null(prefix)) {
  prefix <- if (length(deg_files) == 1) {
    sub("\\.csv$", "", basename(deg_files[1]))
  } else {
    "multi_DEG"
  }
}

# Support comma-separated list of CPM files
count_files <- trimws(strsplit(opt$norm_counts, ",")[[1]])
cat(sprintf("Counts file(s): %d file(s)\n", length(count_files)))

meta_cols <- c("gene_id", "gene_name", "gene_type")

counts_list <- lapply(count_files, function(f) {
  cat(" Reading:", f, "\n")
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(meta_cols, colnames(df))
  if (length(missing) > 0) stop(sprintf("File %s is missing columns: %s", f, paste(missing, collapse = ", ")))
  df
})

# Merge all files by gene_id (outer join — missing values become NA)
if (length(counts_list) == 1) {
  counts <- counts_list[[1]]
} else {
  # Start from the first file's gene annotation columns
  counts <- counts_list[[1]][, meta_cols]
  for (df in counts_list) {
    sample_cols_i <- setdiff(colnames(df), meta_cols)
    dupe_samples  <- intersect(sample_cols_i, colnames(counts))
    if (length(dupe_samples) > 0) {
      warning(sprintf("Duplicate sample names across files — renaming: %s",
                      paste(dupe_samples, collapse = ", ")))
      # append file index suffix to duplicates in df
      fi <- match(list(df), counts_list)
      colnames(df)[colnames(df) %in% dupe_samples] <-
        paste0(dupe_samples, "_f", fi)
      sample_cols_i <- setdiff(colnames(df), meta_cols)
    }
    counts <- merge(counts, df[, c("gene_id", sample_cols_i)],
                    by = "gene_id", all = TRUE)
  }
  # Re-attach gene_name / gene_type from the first file for any gene not merged
  ann <- do.call(rbind, lapply(counts_list, function(df) df[, meta_cols]))
  ann <- ann[!duplicated(ann$gene_id), ]
  counts$gene_name <- ann$gene_name[match(counts$gene_id, ann$gene_id)]
  counts$gene_type <- ann$gene_type[match(counts$gene_id, ann$gene_id)]
  counts <- counts[, c(meta_cols, setdiff(colnames(counts), meta_cols))]
  cat(sprintf("Merged counts matrix: %d genes x %d samples\n",
              nrow(counts), ncol(counts) - length(meta_cols)))
}

# =============================================================================
# Load and filter DEG files (same logic as functional_analysis.r)
# =============================================================================

cat(sprintf("DEG file(s): %d\n", length(deg_files)))

deg_list <- lapply(deg_files, function(f) {
  cat(" Reading:", f, "\n")
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("  %d genes total\n", nrow(df)))
  df
})

filter_sig <- function(df) {
  if (use_logfc) {
    df[!is.na(df$FDR) & df$FDR < opt$fdr & abs(df$logFC) > opt$logfc, ]
  } else {
    df[!is.na(df$FDR) & df$FDR < opt$fdr, ]
  }
}

sig_list <- lapply(seq_along(deg_files), function(i) {
  s <- filter_sig(deg_list[[i]])
  cat(sprintf("  Significant in %s: %d\n", basename(deg_files[i]), nrow(s)))
  s
})

# Union across files; for genes in multiple comparisons keep the row with max |logFC|
sig_all <- do.call(rbind, sig_list)

if (nrow(sig_all) == 0) {
  cat("No significant DEGs found in any file — no heatmap generated.\n")
  quit(save = "no", status = 0)
}

sig_all <- sig_all[order(abs(sig_all$logFC), decreasing = TRUE), ]
sig     <- sig_all[!duplicated(sig_all$gene_id), ]

threshold_label <- if (use_logfc) {
  sprintf("FDR < %.3g AND |logFC| > %.2g", opt$fdr, opt$logfc)
} else {
  sprintf("FDR < %.3g", opt$fdr)
}

cat(sprintf("Union of significant DEGs (%s): %d\n", threshold_label, nrow(sig)))

# Cap at max_genes by max |logFC| across comparisons
if (nrow(sig) > opt$max_genes) {
  cat(sprintf("Capping at top %d genes by max |logFC|\n", opt$max_genes))
  sig <- head(sig, opt$max_genes)
}

cat("Genes in heatmap:", nrow(sig), "\n")

# =============================================================================
# Build CPM matrix for selected genes
# =============================================================================

sample_cols <- setdiff(colnames(counts), meta_cols)

# Subset to sig gene IDs
mat <- counts[counts$gene_id %in% sig$gene_id, ]

if (nrow(mat) == 0) {
  # Fall back to gene_name matching if gene_id not present in counts
  mat <- counts[counts$gene_name %in% sig$gene_name, ]
}

if (nrow(mat) == 0) {
  cat("WARNING: None of the significant genes found in normalized counts file.\n")
  quit(save = "no", status = 0)
}

# Subset samples if requested
if (!is.null(opt$samples)) {
  selected_samples <- trimws(strsplit(opt$samples, ",")[[1]])
  missing <- setdiff(selected_samples, sample_cols)
  if (length(missing) > 0) {
    cat("WARNING: These samples not found in counts file:", paste(missing, collapse = ", "), "\n")
  }
  sample_cols <- intersect(selected_samples, sample_cols)
}

# Compute gap positions (column indices after which to draw a gap)
gaps_col <- NULL
cluster_cols <- TRUE
if (!is.null(opt$gap_samples)) {
  gap_samples <- trimws(strsplit(opt$gap_samples, ",")[[1]])
  missing_gaps <- setdiff(gap_samples, sample_cols)
  if (length(missing_gaps) > 0) {
    cat("WARNING: --gap_samples: these samples not found in column set:",
        paste(missing_gaps, collapse = ", "), "\n")
  }
  gaps_col <- which(sample_cols %in% gap_samples)
  if (length(gaps_col) == 0) {
    cat("WARNING: --gap_samples produced no valid gap positions; ignoring.\n")
    gaps_col <- NULL
  } else {
    gaps_col <- sort(gaps_col)
    cluster_cols <- FALSE   # gaps only make visual sense with a fixed column order
    cat("Gap after column(s):", paste(gaps_col, collapse = ", "),
        "(", paste(sample_cols[gaps_col], collapse = ", "), ")\n")
  }
}

# Build numeric matrix (log2 CPM + pseudocount)
cpm_mat <- as.matrix(mat[, sample_cols])
rownames(cpm_mat) <- mat$gene_name

log_mat <- log2(cpm_mat + 1)

# Z-score by row
z_mat <- t(scale(t(log_mat)))
# Clip extreme values to avoid colour saturation
z_mat <- pmax(pmin(z_mat, 3), -3)

# Restore row order by logFC (up-regulated at top)
gene_order <- sig$gene_name[sig$gene_name %in% rownames(z_mat)]
z_mat <- z_mat[gene_order, , drop = FALSE]

# =============================================================================
# Column annotation from metadata
# =============================================================================

annotation_col <- NA
if (!is.null(opt$metadata)) {
  meta <- read.table(opt$metadata, header = TRUE, stringsAsFactors = FALSE)
  # Keep only columns present in sample_cols
  meta_sub <- meta[meta$SampleId %in% sample_cols, , drop = FALSE]
  rownames(meta_sub) <- meta_sub$SampleId
  meta_sub <- meta_sub[sample_cols, , drop = FALSE]
  # Use all columns except SampleId for annotation
  ann_cols <- setdiff(colnames(meta_sub), "SampleId")
  annotation_col <- meta_sub[, ann_cols, drop = FALSE]
}

# =============================================================================
# Colour palette
# =============================================================================

heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Dynamic figure height: at least 4 inches, add ~0.15 per gene
n_genes   <- nrow(z_mat)
n_samples <- ncol(z_mat)
fig_height <- max(4, min(0.18 * n_genes + 2, 40))
fig_width  <- max(5, min(0.4  * n_samples + 3, 20))

main_title <- sprintf("%s\n(%s,  n = %d genes)", prefix, threshold_label, n_genes)

# Font size: shrink for large gene counts
fontsize_row <- if (n_genes > 80) 4 else if (n_genes > 40) 6 else 8

# =============================================================================
# Plot
# =============================================================================

pdf_path <- file.path(opt$output, paste0("heatmap_", prefix, ".pdf"))
png_path <- file.path(opt$output, paste0("heatmap_", prefix, ".png"))

save_heatmap <- function(file, width, height) {
  pheatmap(
    z_mat,
    color            = heatmap_colors,
    annotation_col   = annotation_col,
    cluster_rows     = TRUE,
    cluster_cols     = cluster_cols,
    gaps_col         = gaps_col,
    show_rownames    = TRUE,
    show_colnames    = TRUE,
    fontsize_row     = fontsize_row,
    fontsize_col     = 8,
    main             = main_title,
    border_color     = NA,
    filename         = file,
    width            = width,
    height           = height
  )
}

save_heatmap(pdf_path, width = fig_width, height = fig_height)
save_heatmap(png_path, width = fig_width, height = fig_height)

cat("Saved:", pdf_path, "\n")
cat("Saved:", png_path, "\n")
