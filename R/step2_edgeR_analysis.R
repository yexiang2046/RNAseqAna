#' do edgeR analysis
#'
#' @param counts A count matrix
#' @param meta  A dataframe containing sample info wiht at least sampleID and treatments
#' @param sample_column column name for sampleID column
#' @param groups A the columns in meta dataframe with treatments
#' @param logFC Threshold for DEGs
#' @param padj  Adjusted p value for difference
#' @return a list of tables contain comparison results


step2_edgeR_analysis <- function(counts = step1_output, meta = NULL,
                                 sample_column = NULL,
                                 groups = "Treatment",
                                 edgeR_file = "edgeR_object.rds",
                                 cpm_file = "cpm.csv",
                                 contrast_list = con,
                                 logFC = 2, padj = 0.05){
  y <- DGEList(counts, samples = meta, group = meta[, groups])

  message("Please make sure the order matches in count matrix and sample info")
  message("counts columns")
  print(colnames(counts))
  message("sample order in metafile")
  print(meta[, sample_column])
  # filtering
  keep <- filterByExpr(y, group = meta[, groups])
  y <- y[keep, ]
  group <- factor(meta[, groups])
  design <- model.matrix(~0+group)
  y <- estimateDisp(y, design)
  saveRDS(y, file = edgeR_file)
  #  plotMDS(y, labels = c("G", "G", "9", "9", "5", "5",
  #                                    "GFP_CLIP", "GFP_CLIP", "K9_CLIP", "K9_CLIP", "ORF52_CLIP", "ORF52_CLIP",
  #                                    "ORF52_input", "ORF52_input"), top = 5000, pch = 19, cex = 0.5, col = "#2E9FDF")
  fit <- glmQLFit(y, design)
  cpm(y) %>% write.csv(cpm_file)

  con <- contrast_list
  #con <- list(con1)
  qlf <- lapply(con, function(x, fit){glmQLFTest(fit, contrast = x)}, fit = fit)
  qlf <- lapply(qlf, function(x){
    x$table$padj <- p.adjust(x$table$PValue, method = "BH")
    x$table$threshold <- as.factor(ifelse(x$table$padj < padj & abs(x$table$logFC) > logFC,
                                          ifelse(x$table$logFC > logFC, 'Up', 'Down'), 'Not'))
    return(x)})

  return(qlf)
}
