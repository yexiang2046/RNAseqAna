#!/usr/bin/env Rscript

# Functional analysis script for differential expression results
# Performs GO, KEGG, Reactome, and GSEA analyses

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
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
# Helper Function: Remove Overlapping Gene Sets
# =============================================================================

#' Remove highly overlapping gene sets based on Jaccard similarity
#'
#' @param enrich_result enrichResult object from clusterProfiler
#' @param overlap_cutoff Jaccard similarity threshold (default 0.7)
#' @return Filtered enrichResult object
remove_overlapping_genesets <- function(enrich_result, overlap_cutoff = 0.7) {
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    return(enrich_result)
  }

  result_df <- enrich_result@result

  # Extract gene lists for each term
  gene_lists <- strsplit(result_df$geneID, "/")
  names(gene_lists) <- result_df$ID

  # Calculate pairwise Jaccard similarity
  n_terms <- length(gene_lists)
  if (n_terms < 2) {
    return(enrich_result)
  }

  to_remove <- c()

  for (i in 1:(n_terms - 1)) {
    if (names(gene_lists)[i] %in% to_remove) next

    for (j in (i + 1):n_terms) {
      if (names(gene_lists)[j] %in% to_remove) next

      genes_i <- gene_lists[[i]]
      genes_j <- gene_lists[[j]]

      # Calculate Jaccard similarity
      intersection <- length(intersect(genes_i, genes_j))
      union <- length(union(genes_i, genes_j))
      jaccard <- intersection / union

      if (jaccard > overlap_cutoff) {
        # Remove the one with higher p-value (less significant)
        pval_i <- result_df$pvalue[result_df$ID == names(gene_lists)[i]]
        pval_j <- result_df$pvalue[result_df$ID == names(gene_lists)[j]]

        if (pval_i <= pval_j) {
          to_remove <- c(to_remove, names(gene_lists)[j])
        } else {
          to_remove <- c(to_remove, names(gene_lists)[i])
          break  # Move to next i since current i is removed
        }
      }
    }
  }

  # Filter result
  if (length(to_remove) > 0) {
    cat("  Removed", length(to_remove), "overlapping terms (>",
        overlap_cutoff * 100, "% similarity)\n")
    enrich_result@result <- result_df[!result_df$ID %in% to_remove, ]
  }

  return(enrich_result)
}

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
    cat("GO BP enrichment: Found", nrow(ego_bp@result), "enriched terms\n")

    # Remove overlapping gene sets
    ego_bp_filtered <- remove_overlapping_genesets(ego_bp, overlap_cutoff = 0.7)

    # Save filtered results
    write.csv(ego_bp_filtered@result,
              file = file.path(output_dir, paste0("GO_BP_", comparison_name, ".csv")),
              row.names = FALSE)

    # Plot filtered results
    if (nrow(ego_bp_filtered@result) > 0) {
      pdf(file.path(output_dir, paste0("GO_BP_", comparison_name, ".pdf")), width = 10, height = 8)
      print(dotplot(ego_bp_filtered, showCategory = 20, title = paste("GO BP -", comparison_name)))
      dev.off()
    }
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
    cat("GO MF enrichment: Found", nrow(ego_mf@result), "enriched terms\n")

    # Remove overlapping gene sets
    ego_mf_filtered <- remove_overlapping_genesets(ego_mf, overlap_cutoff = 0.7)

    # Save filtered results
    write.csv(ego_mf_filtered@result,
              file = file.path(output_dir, paste0("GO_MF_", comparison_name, ".csv")),
              row.names = FALSE)

    # Plot filtered results
    if (nrow(ego_mf_filtered@result) > 0) {
      pdf(file.path(output_dir, paste0("GO_MF_", comparison_name, ".pdf")), width = 10, height = 8)
      print(dotplot(ego_mf_filtered, showCategory = 20, title = paste("GO MF -", comparison_name)))
      dev.off()
    }
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
    cat("GO CC enrichment: Found", nrow(ego_cc@result), "enriched terms\n")

    # Remove overlapping gene sets
    ego_cc_filtered <- remove_overlapping_genesets(ego_cc, overlap_cutoff = 0.7)

    # Save filtered results
    write.csv(ego_cc_filtered@result,
              file = file.path(output_dir, paste0("GO_CC_", comparison_name, ".csv")),
              row.names = FALSE)

    # Plot filtered results
    if (nrow(ego_cc_filtered@result) > 0) {
      pdf(file.path(output_dir, paste0("GO_CC_", comparison_name, ".pdf")), width = 10, height = 8)
      print(dotplot(ego_cc_filtered, showCategory = 20, title = paste("GO CC -", comparison_name)))
      dev.off()
    }
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

    cat("KEGG enrichment: Found", nrow(kegg@result), "enriched pathways\n")

    # Remove overlapping gene sets
    kegg_filtered <- remove_overlapping_genesets(kegg, overlap_cutoff = 0.7)

    # Save filtered results
    write.csv(kegg_filtered@result,
              file = file.path(output_dir, paste0("KEGG_", comparison_name, ".csv")),
              row.names = FALSE)

    # Plot filtered results
    if (nrow(kegg_filtered@result) > 0) {
      pdf(file.path(output_dir, paste0("KEGG_", comparison_name, ".pdf")), width = 10, height = 8)
      print(dotplot(kegg_filtered, showCategory = 20, title = paste("KEGG Pathways -", comparison_name)))
      dev.off()
    }
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
    cat("Reactome enrichment: Found", nrow(reactome@result), "enriched pathways\n")

    # Remove overlapping gene sets
    reactome_filtered <- remove_overlapping_genesets(reactome, overlap_cutoff = 0.7)

    # Save filtered results
    write.csv(reactome_filtered@result,
              file = file.path(output_dir, paste0("Reactome_", comparison_name, ".csv")),
              row.names = FALSE)

    # Plot filtered results
    if (nrow(reactome_filtered@result) > 0) {
      pdf(file.path(output_dir, paste0("Reactome_", comparison_name, ".pdf")), width = 10, height = 8)
      print(dotplot(reactome_filtered, showCategory = 20, title = paste("Reactome Pathways -", comparison_name)))
      dev.off()
    }
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
  
  # Remove duplicate gene names (keep first occurrence with highest |logFC|)
  dup_idx <- duplicated(names(gene_list), fromLast = FALSE)
  if (sum(dup_idx) > 0) {
    cat("  Removed", sum(dup_idx), "duplicate gene names\n")
    gene_list <- gene_list[!dup_idx]
  }

  # Get HALLMARK gene sets
  hallmark <- msigdbr(species = species_msigdb, category = "H")
  hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
  # Convert to list of character vectors
  hallmark_sets <- lapply(hallmark_sets, function(x) as.character(x))

  # Fix for BiocParallel issue
  library(BiocParallel)
  register(SerialParam())

  # Run fgseaMultilevel
  fgsea_results <- fgseaMultilevel(pathways = hallmark_sets,
                                   stats = gene_list,
                                   minSize = 15,
                                   maxSize = 500)

  # Extract key columns and order by padj, excluding list-type columns
  fgsea_output <- data.frame(
    pathway = fgsea_results$pathway,
    pval = fgsea_results$pval,
    padj = fgsea_results$padj,
    log2err = fgsea_results$log2err,
    ES = fgsea_results$ES,
    NES = fgsea_results$NES,
    size = fgsea_results$size,
    stringsAsFactors = FALSE
  )
  fgsea_output <- fgsea_output[order(fgsea_output$padj), ]

  if (nrow(fgsea_output) > 0) {
    write.csv(fgsea_output,
              file = file.path(output_dir, paste0("GSEA_HALLMARK_", comparison_name, ".csv")),
              row.names = FALSE)

    # Plot top 20 pathways
    top_pathways <- head(fgsea_output, 20)

    pdf(file.path(output_dir, paste0("GSEA_HALLMARK_", comparison_name, ".pdf")),
        width = 12, height = 10)

    # Create barplot
    plot_data <- data.frame(
      pathway = gsub("HALLMARK_", "", top_pathways$pathway),
      NES = top_pathways$NES,
      padj = top_pathways$padj,
      stringsAsFactors = FALSE
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

    cat("GSEA HALLMARK: Analyzed", nrow(fgsea_output), "gene sets\n")
    cat("GSEA HALLMARK: Found", sum(fgsea_output$padj < 0.05, na.rm = TRUE), "significant gene sets\n")
  } else {
    cat("GSEA HALLMARK: No results generated\n")
  }
}, error = function(e) {
  cat("Error in GSEA HALLMARK:", e$message, "\n")
})

# =============================================================================
# 5. Over-Representation Analysis (ORA) with HALLMARK gene sets
# =============================================================================

cat("\n=== Running ORA with HALLMARK gene sets ===\n")

tryCatch({
  # Get HALLMARK gene sets
  hallmark <- msigdbr(species = species_msigdb, category = "H")

  # Prepare term2gene format for enricher()
  hallmark_t2g <- hallmark[, c("gs_name", "gene_symbol")]
  colnames(hallmark_t2g) <- c("term", "gene")

  # Use gene symbols from significant DEGs
  ora_hallmark <- enricher(gene = gene_symbols,
                          TERM2GENE = hallmark_t2g,
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

  if (!is.null(ora_hallmark) && nrow(ora_hallmark@result) > 0) {
    cat("ORA HALLMARK: Found", nrow(ora_hallmark@result), "enriched gene sets\n")

    # Remove overlapping gene sets
    ora_hallmark_filtered <- remove_overlapping_genesets(ora_hallmark, overlap_cutoff = 0.7)

    # Save filtered results
    write.csv(ora_hallmark_filtered@result,
              file = file.path(output_dir, paste0("ORA_HALLMARK_", comparison_name, ".csv")),
              row.names = FALSE)

    # Plot filtered results
    if (nrow(ora_hallmark_filtered@result) > 0) {
      pdf(file.path(output_dir, paste0("ORA_HALLMARK_", comparison_name, ".pdf")),
          width = 10, height = 8)
      print(dotplot(ora_hallmark_filtered, showCategory = 20,
                   title = paste("ORA HALLMARK -", comparison_name)))
      dev.off()
    }
  } else {
    cat("ORA HALLMARK: No significant gene sets found\n")
  }
}, error = function(e) {
  cat("Error in ORA HALLMARK:", e$message, "\n")
})

cat("\nFunctional analysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")