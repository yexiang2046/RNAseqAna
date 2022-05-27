
ensembl2gene <- function(ensembl_list){
  # to be add
}



plot_MDS <- function(edger_obj = "edgeR_object.rds", labels = NULL){

  y <- readRDS(edger_obj)
  if(is.null(labels)){
    lab = colnames(y$counts)}
  else {
    lab = labels}
  plotMDS(y, top = 2000, labels = lab, cex = 0.6)
  #text(c(-2.8, -2.8, -2.6, -2.70, -2.3, -2.5, 5.0, 4, 1.2, 1.6, 1.0, 0.5, -1.2, -1.2),
  #     c(0.1, 0.6, 0.1, 0.7, 0.3, 0.6, -1.8, 5.7, -1.9, -1.4, -0.9, -1.2, -0.8, -0.4),
  #     c("G", "G", "9", "9", "5", "5",
  #       "GFP_CLIP", "GFP_CLIP", "K9_CLIP",
  #       "K9_CLIP", "ORF52_CLIP", "ORF52_CLIP",
  #       "ORF52_input", "ORF52_input"),
  #     cex=0.6, pos=4, col="red")
}




plot_vol <- function(qlf, log_FC){
  dat <- qlf$table
  dat$id <- rownames(dat)
  p <- ggplot(data = dat, aes(x = logFC,
                              y = -log10(PValue), label = id,
                              colour = threshold))+
    geom_point() +
    xlab("log2 fold change") +
    ylab("-log10 p-value") +
    ggtitle(gsub("Sample.Description", "", qlf$comparison)) +
    scale_color_manual(values = c("blue", "black", "red")) +
    theme_cowplot() +
    geom_text_repel(data = subset(dat, padj < 0.002 & logFC > log_FC))
  return(p)
}


step5_export_KEGG_terms <- function(edgeR_tb = step4_output, gene_col = 20, comparison = "Name"){
  if("padj" %in% colnames(edgeR_tb)){
    gene_list <- edgeR_tb[edgeR_tb$padj < 0.05,]
    gene_list <- gene_list[!is.na(gene_list$entrez),]
  }else {
    gene_list <- edgeR_tb[!is.na(edgeR_tb$entrez), ]
  }

  kegg_gene_sets = msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
  #m2tg <- kegg_gene_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

  #ekg <- enricher(gene = gene_list$entrez, TERM2GENE = m2tg)

  ekg <- enrichGO(gene = gene_list$entrez, OrgDb = org.Hs.eg.db, ont = "BP")
  ekg <- simplify(ekg, cutoff = 0.7)
  #dotplot(ekg)
  gene_list
  return(ekg)
}


