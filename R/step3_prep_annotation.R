#' prep annotation info
#' @param sp species to retrieve gene annotations
#' @param info_types columns to get, like entrez, symbol
#' @param IdType keys to retrieve info
#' @return A table containing annotation for genes



step3_prep_annotation <- function(sp = "human", IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
                                  info_types = c("entrez", "symbol", "ensembl"),
                                  IdType = "ensembl"){
  src <- src_ucsc(organism = sp)
  IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
  sym <- select(src, keys = IdList, columns = info_types,
                keytype = IdType)
  colnames(sym) <- c("ensembl", "entrez", "symbol")
  return(sym)
}
