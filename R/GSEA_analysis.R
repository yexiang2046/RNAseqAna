#' Perform GSEA analysis



GSEA_enrich <- function(ranked_genes, msigdbr_list){
  fgseaRes <- fgsea(pathways = msigdbr_list, ranked_genes, nPermSimple = 1000000)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
                gseaParam=0.5)
  return(fgseaRes)
}



