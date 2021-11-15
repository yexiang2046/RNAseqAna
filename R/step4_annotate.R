#' Annotate the edgeR comparison tables with more gene info
#' @param qlf A edgeR comparison result
#' @param annotation  Annotation prepared
#' @param Id_col The column used to merge tables
#' @param file_name The filename prefix for edgeR result with annotations
#' @return A table containing edgeR result and annotations





step4_annotate <- function(qlf = step2_output[[1]], annotation = step3_prep_annotation_output, Id_col = "ensembl", file_name = "Name"){
  tb <- qlf$table %>% rownames_to_column()
  tb$rowname <- sapply(tb$rowname, function(x){strsplit(x, "\\.")}[[1]][[1]])
  colnames(tb)[1] <- Id_col
  tb <- left_join(tb, annotation, by = Id_col)
  tb <- tb[!duplicated(tb$ensembl),]
  tb %>% arrange(padj) %>% write_csv(file = paste0(file_name, sep = "_", "DEGs.csv"))
  return(tb)
}
