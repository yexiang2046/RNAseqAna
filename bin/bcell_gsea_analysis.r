#!/usr/bin/env Rscript

# =============================================================================
# B-cell focused GSEA analysis script for differential expression results
# Performs GSEA using:
#   - MSigDB HALLMARK (background)
#   - C2 Canonical Pathways (FULL)
#   - C2 B-cell related pathways (FILTERED)
#   - C7 Immunologic Signatures (B-cell relevant)
#
# Updated: Handles super-long pathway names (Reactome/KEGG/C2) in plots
#          using stringr::str_wrap + taller plots + better y-axis formatting
# =============================================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(BiocParallel)
  library(stringr)          # ← NEW: for wrapping long pathway names
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: bcell_gsea_analysis.r <deg_file> <output_dir> <species> <comparison_name> [use_logfc_filter] [sig_only]\n",
       "       use_logfc_filter: TRUE/1/yes  → FDR < 0.05 AND |logFC| > 1 (DEFAULT)\n",
       "                         FALSE/0/no → FDR < 0.05 only\n",
       "       sig_only:         TRUE/1/yes → save only significant results (padj < 0.05)\n",
       "                         FALSE/0/no → save all results (DEFAULT)")
}

deg_file <- args[1]
output_dir <- args[2]
species <- args[3]
comparison_name <- args[4]

if (length(args) >= 5) {
  use_logfc <- !tolower(args[5]) %in% c("false", "0", "no", "n")
} else {
  use_logfc <- TRUE
}

if (length(args) >= 6) {
  sig_only <- tolower(args[6]) %in% c("true", "1", "yes", "y")
} else {
  sig_only <- FALSE
}

# Set species-specific parameters
if (tolower(species) == "human") {
  species_msigdb <- "Homo sapiens"
} else if (tolower(species) == "mouse") {
  species_msigdb <- "Mus musculus"
} else {
  stop(paste("Unsupported species:", species, ". Use 'human' or 'mouse'."))
}

cat("Starting B-cell focused GSEA for:", comparison_name, "\n")
cat("Species:", species, "\n")
cat("logFC filter enabled (FDR < 0.05 AND |logFC| > 1):", use_logfc, "\n")
cat("Significant results only (padj < 0.05):", sig_only, "\n")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load differential expression results
de_results <- read.csv(deg_file, stringsAsFactors = FALSE)

# =============================================================================
# Prepare ranked gene list for GSEA (uses ALL genes, sorted by logFC)
# =============================================================================
cat("\nPreparing ranked gene list for GSEA...\n")

de_results_complete <- de_results[!is.na(de_results$logFC) &
                                  !is.na(de_results$gene_name) &
                                  de_results$gene_name != "", ]

gene_list <- de_results_complete$logFC
names(gene_list) <- de_results_complete$gene_name
gene_list <- sort(gene_list, decreasing = TRUE)

# Remove duplicate gene names (keep first occurrence)
dup_idx <- duplicated(names(gene_list), fromLast = FALSE)
if (sum(dup_idx) > 0) {
  cat("  Removed", sum(dup_idx), "duplicate gene names\n")
  gene_list <- gene_list[!dup_idx]
}

cat("Ranked gene list ready:", length(gene_list), "genes\n")

# Fix for parallel processing
register(SerialParam())

# =============================================================================
# Helper: Run GSEA + save results + plot top pathways
# =============================================================================
run_bcell_gsea <- function(pathways_list, collection_name, output_prefix, minSize = 15, maxSize = 500) {
  cat("\n=== Running GSEA with", collection_name, "===\n")
  
  fgsea_results <- fgseaMultilevel(pathways = pathways_list,
                                   stats = gene_list,
                                   minSize = minSize,
                                   maxSize = maxSize,
                                   nPermSimple = 10000)

  if (nrow(fgsea_results) == 0) {
    cat(collection_name, ": No results generated\n")
    return(NULL)
  }

  # Clean output
  fgsea_output <- data.frame(
    pathway = fgsea_results$pathway,
    pval = fgsea_results$pval,
    padj = fgsea_results$padj,
    log2err = fgsea_results$log2err,
    ES = fgsea_results$ES,
    NES = fgsea_results$NES,
    size = fgsea_results$size,
    leadingEdge = sapply(fgsea_results$leadingEdge, paste, collapse = ";"),
    stringsAsFactors = FALSE
  )
  fgsea_output <- fgsea_output[order(fgsea_output$padj), ]

  # Apply sig_only filter before saving and plotting
  fgsea_save <- if (sig_only) {
    fgsea_output[!is.na(fgsea_output$padj) & fgsea_output$padj < 0.05, ]
  } else {
    fgsea_output
  }

  write.csv(fgsea_save,
            file = file.path(output_dir, paste0(output_prefix, "_", comparison_name, ".csv")),
            row.names = FALSE)

  n_sig <- sum(fgsea_output$padj < 0.05, na.rm = TRUE)
  cat(collection_name, ": Analyzed", length(pathways_list), "gene sets\n")
  cat(collection_name, ": Found", n_sig, "significant gene sets (padj < 0.05)\n")
  cat(collection_name, ": Saved", nrow(fgsea_save), "gene sets",
      if (sig_only) "(padj < 0.05 only)" else "(all)", "\n")

  # =============================================================================
  # Plot top 20 pathways – now handles SUPER LONG names
  # =============================================================================
  if (nrow(fgsea_save) > 0) {
    top_pathways <- head(fgsea_save[order(abs(fgsea_save$NES), decreasing = TRUE), ], 20)
    n_pathways <- nrow(top_pathways)

    # Scale dimensions to the number of bars
    plot_height <- max(4, 0.45 * n_pathways + 2)
    plot_width  <- if (n_pathways <= 5) 10 else if (n_pathways <= 10) 12 else 14
    # Shrink y-axis text for dense plots
    y_text_size <- if (n_pathways <= 5) 11 else if (n_pathways <= 10) 9.5 else 8.5

    pdf(file.path(output_dir, paste0(output_prefix, "_", comparison_name, ".pdf")),
        width = plot_width, height = plot_height)

    plot_data <- data.frame(
      pathway = gsub("HALLMARK_|REACTOME_|KEGG_|BIOCARTA_", "", top_pathways$pathway),
      NES = top_pathways$NES,
      padj = top_pathways$padj,
      stringsAsFactors = FALSE
    )

    plot_data$pathway <- str_wrap(plot_data$pathway, width = 45)
    plot_data$pathway <- factor(plot_data$pathway, levels = rev(plot_data$pathway))

    p <- ggplot(plot_data, aes(x = NES, y = pathway, fill = padj)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue", limits = c(0, 0.05)) +
      labs(title = paste(collection_name, "-", comparison_name),
           x = "Normalized Enrichment Score (NES)",
           y = "Pathway") +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.y = element_text(size = y_text_size, lineheight = 0.9),
        plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "mm"),
        panel.grid.major.y = element_line(color = "grey90")
      )
    print(p)
    dev.off()
    cat(collection_name, ": Plot saved (", n_pathways, "pathways,",
        round(plot_width, 1), "x", round(plot_height, 1), "in)\n")
  }
  
  return(fgsea_output)
}

# =============================================================================
# 1. GSEA with HALLMARK (background)
# =============================================================================
hallmark <- msigdbr(species = species_msigdb, collection = "H")
hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
hallmark_sets <- lapply(hallmark_sets, as.character)

run_bcell_gsea(hallmark_sets, "HALLMARK", "GSEA_HALLMARK")

# =============================================================================
# 2. GSEA with C2 Canonical Pathways (FULL)
# =============================================================================
c2 <- msigdbr(species = species_msigdb, category = "C2")
c2_sets <- split(c2$gene_symbol, c2$gs_name)
c2_sets <- lapply(c2_sets, as.character)

run_bcell_gsea(c2_sets, "C2_CANONICAL", "GSEA_C2")

# =============================================================================
# 3. GSEA with C2 B-cell related pathways (FILTERED)
# =============================================================================
c2_bcell_keywords <- c("B_CELL_RECEPTOR", "BCR", "B_CELL", "NF_KAPPA_B_IN_B_CELLS",
                       "SIGNALING_BY_THE_B_CELL_RECEPTOR", "ANTIGEN_ACTIVATES_B_CELL",
                       "B_LYMPHOCYTE", "IMMUNOGLOBULIN", "PLASMA_CELL", "CLASS_SWITCH",
                       "SOMATIC_HYPERMUTATION", "MHC_CLASS_II", "B_CELL_ACTIVATION")

c2_bcell <- c2[grepl(paste(c2_bcell_keywords, collapse = "|"), c2$gs_name, ignore.case = TRUE), ]

cat("C2 B-cell filtered gene sets:", nrow(c2_bcell), "rows (", length(unique(c2_bcell$gs_name)), "unique pathways)\n")

if (nrow(c2_bcell) > 0) {
  c2_bcell_sets <- split(c2_bcell$gene_symbol, c2_bcell$gs_name)
  c2_bcell_sets <- lapply(c2_bcell_sets, as.character)
  run_bcell_gsea(c2_bcell_sets, "C2_BCELL_RELEVANT", "GSEA_C2_BCELL")
} else {
  cat("WARNING: No B-cell relevant C2 gene sets found for this species.\n")
}

# =============================================================================
# 4. GSEA with C7 Immunologic Signatures (B-cell relevant)
# =============================================================================
c7 <- msigdbr(species = species_msigdb, category = "C7")

bcell_keywords <- c("B_CELL", "BCELL", "B-CELL", "PLASMA_CELL", "MEMORY_B", 
                    "NAIVE_B", "GERMINAL_CENTER", "BCR", "B_CELL_RECEPTOR", 
                    "IMMUNOGLOBULIN", "CLASS_SWITCH", "SOMATIC_HYPERMUTATION",
                    "B1_CELL", "B-1", "MARGINAL_ZONE")

bcell_c7 <- c7[grepl(paste(bcell_keywords, collapse = "|"), c7$gs_name, ignore.case = TRUE), ]

cat("C7 B-cell filtered gene sets:", nrow(bcell_c7), "rows (", length(unique(bcell_c7$gs_name)), "unique pathways)\n")

if (nrow(bcell_c7) > 0) {
  c7_bcell_sets <- split(bcell_c7$gene_symbol, bcell_c7$gs_name)
  c7_bcell_sets <- lapply(c7_bcell_sets, as.character)
  run_bcell_gsea(c7_bcell_sets, "C7_BCELL_RELEVANT", "GSEA_C7_BCELL")
} else {
  cat("WARNING: No B-cell relevant C7 gene sets found for this species.\n")
}

cat("\nB-cell focused GSEA completed successfully!\n")
cat("Results saved to:", output_dir, "\n")
cat("   → GSEA_HALLMARK_*.csv/pdf\n")
cat("   → GSEA_C2_*.csv/pdf\n")
cat("   → GSEA_C2_BCELL_*.csv/pdf\n")
cat("   → GSEA_C7_BCELL_*.csv/pdf\n")
cat("   → CSV/PDF contain", if (sig_only) "significant (padj < 0.05) results only" else "all results", "\n")