#' Read counts from featurecounts output
#'
#' @param featurecount_file A rds file containing featurecounts results
#' @return the count matrix with ensembl ID as rowname

step1_get_counts <- function(featurecount_file = "featurecounts.rds", type = "command_line",
                             col_names = c("Geneid", "Chr", "Start", "End", "Strand", "Length", "GFP1", "GFP2", "LANA1", "LANA2"),
                             counts_col = c("GFP1", "GFP2", "LANA1", "LANA2")){
  if(type == "command_line"){
    counts <- read.table(file = featurecount_file, sep = "\t", header = TRUE, skip = 1)
    colnames(counts) <- col_names
    cts <- counts[, counts_col]
    rownames(cts) <- counts$Geneid
  }else{
    counts <- readRDS(featurecount_file)
    cts <- counts$counts
  }
  return(cts)
}
