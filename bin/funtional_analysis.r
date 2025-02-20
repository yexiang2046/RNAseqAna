#! /usr/bin/env Rscript

# this script performs functional analysis on the differentially expressed genes
# it uses the clusterProfiler package to perform gene set enrichment analysis
# it uses fgsea package to perform gene set enrichment analysis
# the database used is msigdb

# load the necessary libraries
library(clusterProfiler)
library(fgsea)
library(msigdbr)

# load the msigdb database
msigdb <- msigdbr(species = "Homo sapiens", category = "C2")
# prepare genesets for fgsea
genesets <- split(msigdb$gene_symbol, msigdb$gs_name)

# get list of DEG files
deg_files <- list.files(path = "de_results", pattern = "DEG_.*\\.csv", full.names = TRUE)

# iterate over each DEG file
for (deg_file in deg_files) {
  # extract the comparison name from the file name
  comparison_name <- sub(".*DEG_(.*)\\.csv", "\\1", basename(deg_file))
  
  # load the differential expression results
  de_results <- read.csv(deg_file)
  
  # perform gene set enrichment analysis
  ego <- enrichKEGG(gene = de_results$gene,
                    universe = de_results$gene,
                    organism = "hsa",
                    pAdjustMethod = "BH")
  # plot the results
  pdf(paste0("ego_", comparison_name, ".pdf"))
  dotplot(ego, showCategory = 10)
  dev.off()
  
  # perform fgsea analysis
  fgsea_results <- fgsea(pathways = genesets, stats = de_results$log2FoldChange)
  # select the top 10 pathways upregulated and downregulated
  fgsea_results_top10 <- fgsea_results[order(fgsea_results$padj),][1:10,]

  # plot the results
  pdf(paste0("fgsea_", comparison_name, "top10.pdf"))
  plotGseaTable(genesets, fgsea_results_top10, gseaParam = 0.5)
  dev.off()
  
  # save the results
  write.csv(ego, file = paste0("ego_", comparison_name, ".csv"))
  write.csv(fgsea_results, file = paste0("fgsea_", comparison_name, ".csv"))
}