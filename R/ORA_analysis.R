#' Perform ORA analysis
#' @param genelist a list of gene entrez ID or symbol
#' @param TREM2genes a dataframe of terms and gene entrez ID or symbol
#' @return A dataframe of enriched terms for the gene list

ORA_enrich <- function(genelist, TERM2genes){
  gene_list <- genelist

  
  ekg <- enricher(gene = gene_list$entrez, TERM2GENE = m2tg)
  #dotplot(ekg)
  gene_list
  return(ekg)
}