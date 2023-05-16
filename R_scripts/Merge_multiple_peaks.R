library(rtracklayer)
library(GenomicRanges)
library(dplyr)

T48hr_ac_peaks <- c("T48hr_CB839_H3K27ac_idr.05.txt", "T48hr_CTL_H3K27ac_idr.05.txt", "T48hr_DON_H3K27ac_idr.05.txt", "T48hr_NoQ_H3K27ac_idr.05.txt")
T48hr_ac_peaks <- paste0("TXT/", T48hr_ac_peaks)
T48hr_me3_peaks <- c("T48hr_CB839_H3K27me3_idr.05.txt", "T48hr_CTL_H3K27me3_idr.05.txt", "T48hr_DON_H3K27me3_idr.05.txt", "T48hr_NoQ_H3K27me3_idr.05.txt")
T48hr_me3_peaks <- paste0("TXT/", T48hr_me3_peaks)

# functions to turn idr output into genomic ranges
import_idr <- function(x){
    out <- read.delim(file = x, header = FALSE)
    out <- out[,c(1:7,11)]
    # turn the fold enrichment into mean for the replicates
    out$V7 <- out$V7/2
    names(out) <- c("Chr", "Start", "End", "Name", "Score", "Strand", "FoldEnrichment", "-log10(FDR)")
    out_ig <- makeGRangesFromDataFrame(out, keep.extra.columns = TRUE)
}


T48hr_ac <- lapply(T48hr_ac_peaks, function(x){import_idr(x)})

T48hr_ac <- GRangesList(T48hr_ac)
mcols(T48hr_ac) <- c("CTL", "CB839", "DON", "NoQ")


T48hr_ac_gr <- unlist(T48hr_ac)

mcols(T48hr_ac_gr)$Name <- c(paste0("CTL_", 1:length(T48hr_ac[[1]])), paste0("CB839_", 1:length(T48hr_ac[[2]])), paste0("DON_", 1:length(T48hr_ac[[3]])), paste0("NoQ_", 1:length(T48hr_ac[[4]])))


T48hr_ac_reduce <- reduce(T48hr_ac_gr, ignore.strand = TRUE, drop.empty.ranges = FALSE, with.revmap = TRUE)
revmap <- mcols(T48hr_ac_reduce)$revmap
# relist(T48hr_ac_gr[unlist(revmap)], revmap)
##
T48hr_ac_expandgr <- relist(mcols(T48hr_ac_gr)[unlist(revmap), ], revmap)

T48hr_ac_df <- data.frame(T48hr_ac_reduce)

for (i in 1:length(T48hr_ac_expandgr)) {
    T48hr_ac_df[i, c("CTL_FoldEnrichment", "CTL_-log10(FDR)")] <- as.vector(unlist(T48hr_ac_expandgr[[i]][grep("CTL", T48hr_ac_expandgr[[i]]$Name),]))[3:4]
    T48hr_ac_df[i, c("CB839_FoldEnrichment", "CB839_-log10(FDR)")] <- as.vector(unlist(T48hr_ac_expandgr[[i]][grep("CB839", T48hr_ac_expandgr[[i]]$Name),]))[3:4]
    T48hr_ac_df[i, c("DON_FoldEnrichment", "DON_-log10(FDR)")] <- as.vector(unlist(T48hr_ac_expandgr[[i]][grep("DON", T48hr_ac_expandgr[[i]]$Name),]))[3:4]
    T48hr_ac_df[i, c("NoQ_FoldEnrichment", "NoQ_-log10(FDR)")] <- as.vector(unlist(T48hr_ac_expandgr[[i]][grep("NoQ", T48hr_ac_expandgr[[i]]$Name),]))[3:4]
}

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(Organism.dplyr)
library(tidyverse)




T48hr_ac_meta_gr <- makeGRangesFromDataFrame(T48hr_ac_df, keep.extra.columns = TRUE)
# T48hr_ac_meta_gr$score <- T48hr_ac_reduce$score
# export.bed(T48hr_ac_meta_gr, con = "T48hr_ac_conditions_stats.bed")
T48hr_ac_meta_gr_anno <- annotatePeak(T48hr_ac_meta_gr, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)

src <- src_ucsc(organism = "mouse", genome = "mm10")
IdList <- T48hr_ac_meta_gr_anno@anno$geneId

sym <- select(src, keys = IdList, columns = c("entrez", "symbol"),
                keytype = "entrez")

colnames(sym) <- c("geneId", "symbol")
sym <- sym[!duplicated(sym$geneId),]

T48hr_ac_anno <- T48hr_ac_meta_gr_anno@anno
T48hr_ac_anno_df <- data.frame(T48hr_ac_anno)[, c(1:5,7:23)]
T48hr_ac_anno_df <- left_join(T48hr_ac_anno_df, sym, by = "geneId")
T48hr_ac_anno_df[is.na(T48hr_ac_anno_df)] <- "Not Available"
T48hr_ac_anno_df %>% data.frame() %>% write.csv(file = "T48hr_ac_annotation.csv", row.names = FALSE)



# for H3K27me3
T48hr_me3 <- lapply(T48hr_me3_peaks, function(x){import_idr(x)})

T48hr_me3 <- GRangesList(T48hr_me3)
mcols(T48hr_me3) <- c("CTL", "CB839", "DON", "NoQ")


T48hr_me3_gr <- unlist(T48hr_me3)

mcols(T48hr_me3_gr)$Name <- c(paste0("CTL_", 1:length(T48hr_me3[[1]])), paste0("CB839_", 1:length(T48hr_me3[[2]])), paste0("DON_", 1:length(T48hr_me3[[3]])), paste0("NoQ_", 1:length(T48hr_me3[[4]])))


T48hr_me3_reduce <- reduce(T48hr_me3_gr, ignore.strand = TRUE, drop.empty.ranges = FALSE, with.revmap = TRUE)
revmap <- mcols(T48hr_me3_reduce)$revmap
# relist(T48hr_me3_gr[unlist(revmap)], revmap)
##
T48hr_me3_expandgr <- relist(mcols(T48hr_me3_gr)[unlist(revmap), ], revmap)

T48hr_me3_df <- data.frame(T48hr_me3_reduce)

for (i in 1:length(T48hr_me3_expandgr)) {
    T48hr_me3_df[i, c("CTL_FoldEnrichment", "CTL_-log10(FDR)")] <- as.vector(unlist(T48hr_me3_expandgr[[i]][grep("CTL", T48hr_me3_expandgr[[i]]$Name),]))[3:4]
    T48hr_me3_df[i, c("CB839_FoldEnrichment", "CB839_-log10(FDR)")] <- as.vector(unlist(T48hr_me3_expandgr[[i]][grep("CB839", T48hr_me3_expandgr[[i]]$Name),]))[3:4]
    T48hr_me3_df[i, c("DON_FoldEnrichment", "DON_-log10(FDR)")] <- as.vector(unlist(T48hr_me3_expandgr[[i]][grep("DON", T48hr_me3_expandgr[[i]]$Name),]))[3:4]
    T48hr_me3_df[i, c("NoQ_FoldEnrichment", "NoQ_-log10(FDR)")] <- as.vector(unlist(T48hr_me3_expandgr[[i]][grep("NoQ", T48hr_me3_expandgr[[i]]$Name),]))[3:4]
}

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(Organism.dplyr)
library(tidyverse)




T48hr_me3_meta_gr <- makeGRangesFromDataFrame(T48hr_me3_df, keep.extra.columns = TRUE)
# T48hr_me3_meta_gr$score <- T48hr_me3_reduce$score
# export.bed(T48hr_me3_meta_gr, con = "T48hr_me3_conditions_stats.bed")
T48hr_me3_meta_gr_anno <- annotatePeak(T48hr_me3_meta_gr, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)

src <- src_ucsc(organism = "mouse", genome = "mm10")
IdList <- T48hr_me3_meta_gr_anno@anno$geneId

sym <- select(src, keys = IdList, columns = c("entrez", "symbol"),
                keytype = "entrez")

colnames(sym) <- c("geneId", "symbol")
sym <- sym[!duplicated(sym$geneId),]

T48hr_me3_anno <- T48hr_me3_meta_gr_anno@anno
T48hr_me3_anno_df <- data.frame(T48hr_me3_anno)[, c(1:5,7:23)]
T48hr_me3_anno_df <- left_join(T48hr_me3_anno_df, sym, by = "geneId")
T48hr_me3_anno_df[is.na(T48hr_me3_anno_df)] <- "Not Available"
T48hr_me3_anno_df %>% data.frame() %>% write.csv(file = "T48hr_me3_annotation.csv", row.names = FALSE)
