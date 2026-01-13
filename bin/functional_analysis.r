#!/usr/bin/env Rscript

# Functional analysis script for differential expression results
# Performs GO, KEGG, Reactome, and GSEA analyses

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(fgsea)
  library(msigdbr)
  library(enrichplot)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: functional_analysis.r <deg_file> <output_dir> <species> <comparison_name>")
}

deg_file <- args[1]
output_dir <- args[2]
species <- args[3]
comparison_name <- args[4]

# Set species-specific parameters
if (tolower(species) == "human") {
  species_msigdb <- "Homo sapiens"
  organism_kegg <- "hsa"
  orgdb <- org.Hs.eg.db
} else if (tolower(species) == "mouse") {
  species_msigdb <- "Mus musculus"
  organism_kegg <- "mmu"
  orgdb <- org.Mm.eg.db
} else {
  stop(paste("Unsupported species:", species, ". Use 'human' or 'mouse'."))
}

cat("Starting functional analysis for:", comparison_name, "\n")
cat("Species:", species, "\n")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load differential expression results
de_results <- read.csv(deg_file, stringsAsFactors = FALSE)

# Filter significant DEGs (FDR < 0.05, |logFC| > 1)
sig_genes <- de_results[!is.na(de_results$FDR) &
                        de_results$FDR < 0.05 &
                        abs(de_results$logFC) > 1, ]

cat("Total genes:", nrow(de_results), "\n")
cat("Significant DEGs:", nrow(sig_genes), "\n")

# Check if we have enough genes
if (nrow(sig_genes) < 5) {
  cat("WARNING: Less than 5 significant DEGs. Skipping enrichment analysis.\n")
  quit(save = "no", status = 0)
}

# Prepare gene lists
gene_symbols <- sig_genes$gene_name
gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]

if (length(gene_symbols) < 5) {
  cat("WARNING: Less than 5 genes with valid symbols. Skipping enrichment analysis.\n")
  quit(save = "no", status = 0)
}

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(gene_symbols,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = orgdb)

if (nrow(gene_entrez) < 5) {
  cat("WARNING: Less than 5 genes with valid Entrez IDs. Skipping enrichment analysis.\n")
  quit(save = "no", status = 0)
}

cat("Genes mapped to Entrez IDs:", nrow(gene_entrez), "\n")

# =============================================================================
# 1. Gene Ontology (GO) Enrichment Analysis
# =============================================================================

cat("\n=== Running GO Enrichment Analysis ===\n")

# GO Biological Process
tryCatch({
  ego_bp <- enrichGO(gene = gene_entrez$ENTREZID,
                     OrgDb = orgdb,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)

  if (!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
    write.csv(ego_bp@result,
              file = file.path(output_dir, paste0("GO_BP_", comparison_name, ".csv")),
              row.names = FALSE)

    pdf(file.path(output_dir, paste0("GO_BP_", comparison_name, ".pdf")), width = 10, height = 8)
    print(dotplot(ego_bp, showCategory = 20, title = paste("GO BP -", comparison_name)))
    dev.off()
    cat("GO BP enrichment: Found", nrow(ego_bp@result), "enriched terms\n")
  } else {
    cat("GO BP enrichment: No significant terms found\n")
  }
}, error = function(e) {
  cat("Error in GO BP enrichment:", e$message, "\n")
})

# GO Molecular Function
tryCatch({
  ego_mf <- enrichGO(gene = gene_entrez$ENTREZID,
                     OrgDb = orgdb,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)

  if (!is.null(ego_mf) && nrow(ego_mf@result) > 0) {
    write.csv(ego_mf@result,
              file = file.path(output_dir, paste0("GO_MF_", comparison_name, ".csv")),
              row.names = FALSE)

    pdf(file.path(output_dir, paste0("GO_MF_", comparison_name, ".pdf")), width = 10, height = 8)
    print(dotplot(ego_mf, showCategory = 20, title = paste("GO MF -", comparison_name)))
    dev.off()
    cat("GO MF enrichment: Found", nrow(ego_mf@result), "enriched terms\n")
  } else {
    cat("GO MF enrichment: No significant terms found\n")
  }
}, error = function(e) {
  cat("Error in GO MF enrichment:", e$message, "\n")
})

# GO Cellular Component
tryCatch({
  ego_cc <- enrichGO(gene = gene_entrez$ENTREZID,
                     OrgDb = orgdb,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)

  if (!is.null(ego_cc) && nrow(ego_cc@result) > 0) {
    write.csv(ego_cc@result,
              file = file.path(output_dir, paste0("GO_CC_", comparison_name, ".csv")),
              row.names = FALSE)

    pdf(file.path(output_dir, paste0("GO_CC_", comparison_name, ".pdf")), width = 10, height = 8)
    print(dotplot(ego_cc, showCategory = 20, title = paste("GO CC -", comparison_name)))
    dev.off()
    cat("GO CC enrichment: Found", nrow(ego_cc@result), "enriched terms\n")
  } else {
    cat("GO CC enrichment: No significant terms found\n")
  }
}, error = function(e) {
  cat("Error in GO CC enrichment:", e$message, "\n")
})

# =============================================================================
# 2. KEGG Pathway Enrichment
# =============================================================================

cat("\n=== Running KEGG Pathway Enrichment ===\n")

tryCatch({
  kegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                     organism = organism_kegg,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)

  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    # Convert Entrez IDs to gene symbols for readability
    kegg <- setReadable(kegg, OrgDb = orgdb, keyType = "ENTREZID")

    write.csv(kegg@result,
              file = file.path(output_dir, paste0("KEGG_", comparison_name, ".csv")),
              row.names = FALSE)

    pdf(file.path(output_dir, paste0("KEGG_", comparison_name, ".pdf")), width = 10, height = 8)
    print(dotplot(kegg, showCategory = 20, title = paste("KEGG Pathways -", comparison_name)))
    dev.off()
    cat("KEGG enrichment: Found", nrow(kegg@result), "enriched pathways\n")
  } else {
    cat("KEGG enrichment: No significant pathways found\n")
  }
}, error = function(e) {
  cat("Error in KEGG enrichment:", e$message, "\n")
})

# =============================================================================
# 3. Reactome Pathway Enrichment
# =============================================================================

cat("\n=== Running Reactome Pathway Enrichment ===\n")

tryCatch({
  reactome <- enrichPathway(gene = gene_entrez$ENTREZID,
                           organism = ifelse(tolower(species) == "human", "human", "mouse"),
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)

  if (!is.null(reactome) && nrow(reactome@result) > 0) {
    write.csv(reactome@result,
              file = file.path(output_dir, paste0("Reactome_", comparison_name, ".csv")),
              row.names = FALSE)

    pdf(file.path(output_dir, paste0("Reactome_", comparison_name, ".pdf")), width = 10, height = 8)
    print(dotplot(reactome, showCategory = 20, title = paste("Reactome Pathways -", comparison_name)))
    dev.off()
    cat("Reactome enrichment: Found", nrow(reactome@result), "enriched pathways\n")
  } else {
    cat("Reactome enrichment: No significant pathways found\n")
  }
}, error = function(e) {
  cat("Error in Reactome enrichment:", e$message, "\n")
})

# =============================================================================
# 4. Gene Set Enrichment Analysis (GSEA) with HALLMARK gene sets
# =============================================================================

cat("\n=== Running GSEA with HALLMARK gene sets ===\n")

tryCatch({
  # Prepare ranked gene list (use all genes, not just significant ones)
  de_results_complete <- de_results[!is.na(de_results$logFC) &
                                    !is.na(de_results$gene_name) &
                                    de_results$gene_name != "", ]

  # Create ranked list by logFC
  gene_list <- de_results_complete$logFC
  names(gene_list) <- de_results_complete$gene_name
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Get HALLMARK gene sets
  hallmark <- msigdbr(species = species_msigdb, category = "H")
  hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)

  # Run fgsea
  fgsea_results <- fgsea(pathways = hallmark_sets,
                        stats = gene_list,
                        minSize = 15,
                        maxSize = 500,
                        nperm = 10000)

  # Order by padj
  fgsea_results <- fgsea_results[order(fgsea_results$padj), ]

  if (nrow(fgsea_results) > 0) {
    write.csv(fgsea_results,
              file = file.path(output_dir, paste0("GSEA_HALLMARK_", comparison_name, ".csv")),
              row.names = FALSE)

    # Plot top 20 pathways
    top_pathways <- head(fgsea_results, 20)

    pdf(file.path(output_dir, paste0("GSEA_HALLMARK_", comparison_name, ".pdf")),
        width = 12, height = 10)

    # Create barplot
    plot_data <- data.frame(
      pathway = gsub("HALLMARK_", "", top_pathways$pathway),
      NES = top_pathways$NES,
      padj = top_pathways$padj
    )
    plot_data$pathway <- factor(plot_data$pathway, levels = rev(plot_data$pathway))

    p <- ggplot(plot_data, aes(x = NES, y = pathway, fill = padj)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue") +
      labs(title = paste("GSEA HALLMARK -", comparison_name),
           x = "Normalized Enrichment Score (NES)",
           y = "Pathway") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    print(p)
    dev.off()

    cat("GSEA HALLMARK: Analyzed", nrow(fgsea_results), "gene sets\n")
    cat("GSEA HALLMARK: Found", sum(fgsea_results$padj < 0.05, na.rm = TRUE), "significant gene sets\n")
  } else {
    cat("GSEA HALLMARK: No results generated\n")
  }
}, error = function(e) {
  cat("Error in GSEA HALLMARK:", e$message, "\n")
})

cat("\nFunctional analysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")
