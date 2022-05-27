require(dplyr)
require(readr)
require(readxl)
require(edgeR)
require(clusterProfiler)
require(fgsea)
require(msigdbr)
#require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(org.Mm.eg.db)
require(Organism.dplyr)
require(ggfortify) # for PCA plot
require(cowplot)
require(ggrepel)
require(tibble)
require(tidyr)
#require(TCseq)

if(FALSE){
  metaData <- read_xlsx("SAMPLE_SHEET/7342_MM_Samples_IDs.xlsx")
  colnames(metaData) <- c("sampleid", "Description")
  metaData$timepoint <- sapply(metaData$Description, function(x){strsplit(x, split = "\\s")[[1]][1]})
  metaData$treatment <- sapply(metaData$Description, function(x){strsplit(x, split = "\\s")[[1]][2]})
  metaData$treatment[c(1:8, 21:24, 37:40)] <- "Vehicle"
  metaData$treatment[c(17:20,  33:36, 49:52)] <- "No_Q"

  metaData$group <- paste0(metaData$treatment, ".", metaData$timepoint)
  # metaData$group <- c(1,1,1,1,rep(2,16), rep(3,16),rep(4,16))

  bam_files <- c("7342-MM-1.filtered.bam","7342-MM-2.filtered.bam","7342-MM-3.filtered.bam",
                "7342-MM-4.filtered.bam","7342-MM-5.filtered.bam","7342-MM-6.filtered.bam",
                "7342-MM-7.filtered.bam","7342-MM-8.filtered.bam","7342-MM-9.filtered.bam",
                "7342-MM-10.filtered.bam","7342-MM-11.filtered.bam","7342-MM-12.filtered.bam",
                "7342-MM-13.filtered.bam","7342-MM-14.filtered.bam","7342-MM-15.filtered.bam",
                "7342-MM-16.filtered.bam","7342-MM-17.filtered.bam","7342-MM-18.filtered.bam",
                "7342-MM-19.filtered.bam","7342-MM-20.filtered.bam","7342-MM-21.filtered.bam",
                "7342-MM-22.filtered.bam","7342-MM-23.filtered.bam","7342-MM-24.filtered.bam",
                "7342-MM-25.filtered.bam","7342-MM-26.filtered.bam","7342-MM-27.filtered.bam",
                "7342-MM-28.filtered.bam","7342-MM-29.filtered.bam","7342-MM-30.filtered.bam",
                "7342-MM-31.filtered.bam","7342-MM-32.filtered.bam","7342-MM-33.filtered.bam",
                "7342-MM-34.filtered.bam","7342-MM-35.filtered.bam","7342-MM-36.filtered.bam",
                "7342-MM-37.filtered.bam","7342-MM-38.filtered.bam","7342-MM-39.filtered.bam",
                "7342-MM-40.filtered.bam","7342-MM-41.filtered.bam","7342-MM-42.filtered.bam",
                "7342-MM-43.filtered.bam","7342-MM-44.filtered.bam","7342-MM-45.filtered.bam",
                "7342-MM-46.filtered.bam","7342-MM-47.filtered.bam","7342-MM-48.filtered.bam",
                "7342-MM-49.filtered.bam","7342-MM-50.filtered.bam",
                "7342-MM-51.filtered.bam","7342-MM-52.filtered.bam")

  metaData$bamfile <- bam_files

  #CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES <- metaData[,c("sampleid", "timepoint", "group", "bamfile")]
  #CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES <- metaData
  #use_data(CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES, overwrite = TRUE)

  counts <- read.table(file = "counts/CD8_RNAseq_counts.txt", sep = "\t", header = TRUE, skip = 1)
  colnames(counts) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", metaData$sampleid)
  Geneid <- counts$Geneid

  CD8_RNAseq_GLUTAMINE_METABOLISM_FEATURE_COUNT <- counts
  CD8_RNAseq_meta <- metaData

## raw counts
  library(devtools)
  # use_data(CD8_RNAseq_GLUTAMINE_METABOLISM_FEATURE_COUNT, overwrite = TRUE)
  # use_data(CD8_RNAseq_meta)

  counts <- read.table(file = "counts/CD8_RNAseq_counts.txt", sep = "\t", header = TRUE, skip = 1)
  colnames(counts) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", metaData$sampleid)

  src <- src_ucsc(organism = "mouse", genome = "mm10")
  IdList <- sapply(counts$Geneid, function(x){strsplit(x, "\\.")}[[1]][[1]])
  counts$Geneid <- IdList
  sym <- select(src, keys = IdList, columns = c("entrez", "symbol", "ensembl", "gene_chrom", "gene_start", "gene_end"),
                keytype = "ensembl")
  colnames(sym) <- c("Geneid", "entrez", "symbol", "chr", "start", "end")

  counts <- left_join(counts, sym, by = "Geneid")
  counts <- counts[!is.na(counts$entrez), ]
  counts <- counts[!duplicated(counts$entrez), ]
  rownames(counts) <- counts$entrez
  CD8_RNAseq_GLUTAMINE_METABOLISM_COUNTS <- counts

  use_data(CD8_RNAseq_GLUTAMINE_METABOLISM_COUNTS, overwrite = TRUE)
  ctx <- as.matrix(CD8_RNAseq_GLUTAMINE_METABOLISM_COUNTS[, 7:58])
  rownames(ctx) <- CD8_RNAseq_GLUTAMINE_METABOLISM_COUNTS$entrez
  y <- DGEList(ctx, samples = CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES,
              group = CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES$group,
              genes = CD8_RNAseq_GLUTAMINE_METABOLISM_COUNTS[, c("Geneid", "entrez", "symbol", "chr", "start", "end")])

# filtering
  keep <- filterByExpr(y, group = CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES$group)
  y <- y[keep, ]
  y <- calcNormFactors(y)
  group <- factor(metaData$group)
  design <- model.matrix(~0+group)
  y <- estimateDisp(y, design)
  saveRDS(y, file = "RDS/CD8_RNAseq_edgeR_object.rds")

}

data("CD8_RNAseq_meta")
data("CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES")
data("CD8_RNAseq_GLUTAMINE_METABOLISM_COUNTS")



y <- readRDS(file = "RDS/CD8_RNAseq_edgeR_object.rds")
if(FALSE){
  colors <- c(rep("black", 8), rep("green", 4), rep("red", 4), rep("blue", 4),
             rep("black", 4), rep("green", 4), rep("red", 4), rep("blue", 4),
             rep("black", 4), rep("green", 4), rep("red", 4), rep("blue", 4))
  points <- c(rep(0, 8), rep(1, 4), rep(2, 4), rep(23, 4),
             rep(0, 4), rep(1, 4), rep(2, 4), rep(23, 4),
             rep(0, 4), rep(1, 4), rep(2, 4), rep(23, 4))


  plotMDS(y, top = 5000, pch = points, col = colors)
}


if(FALSE){
  library(pheatmap)
  logNormCount <- cpm(y, log = TRUE)
#heatmap(cor(logNormCount))

  hm <- cor(logNormCount)
  ann_col <- data.frame(CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES[, c("timepoint", "treatment")])
  rownames(ann_col) <- CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES$sampleid
  ann_col$timepoint <- factor(ann_col$timepoint, levels = c("0hr", "4hr", "24hr", "48hr"))
  ann_col$treatment <- factor(ann_col$treatment, levels = c("Vehicle", "CB839", "DON", "No_Q"))
  pheatmap(hm, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = ann_col)

  norm_count <- cpm(y) %>% data.frame() %>% rownames_to_column()
  colnames(norm_count)[1] <- "entrez"
  genes <- y$genes
  norm_count <- left_join(norm_count, genes, by = "entrez")
  norm_count <- norm_count[, c(1, 54:58, 2:53)]
  norm_count <- norm_count[!is.na(norm_count$entrez),]
  norm_count %>% write.csv("CD8_RNAseq_cpm.csv")
}

group <- factor(CD8_RNAseq_meta$group)
design <- model.matrix(~0+group)
con1 <- makeContrasts(groupVehicle.4hr - groupVehicle.0hr, levels = design)
con2 <- makeContrasts(groupVehicle.24hr - groupVehicle.0hr, levels = design)
con3 <- makeContrasts(groupVehicle.48hr - groupVehicle.0hr, levels = design)
con4 <- makeContrasts(groupCB839.4hr - groupVehicle.0hr, levels = design)
con5 <- makeContrasts(groupCB839.24hr - groupVehicle.0hr, levels = design)
con6 <- makeContrasts(groupCB839.48hr - groupVehicle.0hr, levels = design)
con7 <- makeContrasts(groupDON.4hr - groupVehicle.0hr, levels = design)
con8 <- makeContrasts(groupDON.24hr - groupVehicle.0hr, levels = design)
con9 <- makeContrasts(groupDON.48hr - groupVehicle.0hr, levels = design)
con10 <- makeContrasts(groupNo_Q.4hr - groupVehicle.0hr, levels = design)
con11 <- makeContrasts(groupNo_Q.24hr - groupVehicle.0hr, levels = design)
con12 <- makeContrasts(groupNo_Q.48hr - groupVehicle.0hr, levels = design)

con13 <- makeContrasts(groupCB839.4hr - groupVehicle.4hr, levels = design)
con14 <- makeContrasts(groupDON.4hr - groupVehicle.4hr, levels = design)
con15 <- makeContrasts(groupNo_Q.4hr - groupVehicle.4hr, levels = design)

con16 <- makeContrasts(groupCB839.24hr - groupVehicle.24hr, levels = design)
con17 <- makeContrasts(groupDON.24hr - groupVehicle.24hr, levels = design)
con18 <- makeContrasts(groupNo_Q.24hr - groupVehicle.24hr, levels = design)

con19 <- makeContrasts(groupCB839.48hr - groupVehicle.48hr, levels = design)
con20 <- makeContrasts(groupDON.48hr - groupVehicle.48hr, levels = design)
con21 <- makeContrasts(groupNo_Q.48hr - groupVehicle.48hr, levels = design)

con <- list(con1, con2, con3, con4, con5, con6, con7, con8, con9, con10, con11,
            con12, con13, con14, con15, con16, con17, con18, con19, con20, con21)
#con <- list(con1)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- lapply(con, function(x, fit){glmQLFTest(fit, contrast = x)}, fit = fit)
qlf <- lapply(qlf, function(x){
  x$table$padj <- p.adjust(x$table$PValue, method = "BH")
  x$table$threshold <- as.factor(ifelse(x$table$padj < 0.05 & abs(x$table$logFC) > 1,
                                        ifelse(x$table$logFC > 1, 'Up', 'Down'), 'Not'))
  return(x)
})

## function
## export DEGs from DGEList
export_tb_from_dgelist <- function(x){
  tb <- x$table %>% data.frame() %>% rownames_to_column()
  #colnames(tb)[1] <- "entrez"
  genes <- x$genes %>% data.frame() %>% rownames_to_column()
  tb <- left_join(tb, genes, by = "rowname")
  file_n <- x$comparison
  file_n <- gsub("\\*group", 'group', file_n)
  file_n <- gsub('-1', 'ref', file_n)
  file_n <- gsub('1', '_for_', file_n)
  file_n <- paste0(file_n, "_DEGs.csv")
  tb[,c(8:13, 2:7)] %>% write_csv(file = file_n)
  return(tb)
}

filter_DEGs <- function(tb){
  tb %>% filter(padj < 0.05)
}

## export resutlt for all DEGs
if(FALSE){
  lapply(qlf, export_tb_from_dgelist)
}


T48hr_CB839_vs_CTL <- export_tb_from_dgelist(qlf[[19]])
T48hr_DON_vs_CTL <- export_tb_from_dgelist(qlf[[20]])
T48hr_NoQ_vs_CTL <- export_tb_from_dgelist(qlf[[21]])

T48hr_CB839_vs_CTL_DEGs <- filter_DEGs(T48hr_CB839_vs_CTL)


## T48hr merged peaks from CUT&RUN-------------------

get_ProPeak_genes <- function(ann_tb){
  anno_tb <- read.delim2(file = ann_tb)
  anno_tb <- anno_tb[abs(anno_tb$Distance.to.TSS) < 2000, ]
  return(anno_tb)
}

T48hr_ac_peaks <- c("T48hr_CTL_H3K27ac_idr.05.annotation.txt", "T48hr_CB839_H3K27ac_idr.05.annotation.txt", "T48hr_DON_H3K27ac_idr.05.annotation.txt", "T48hr_NoQ_H3K27ac_idr.05.annotation.txt")
T48hr_me3_peaks <- c("T48hr_CTL_H3K27me3_idr.05.annotation.txt", "T48hr_CB839_H3K27me3_idr.05.annotation.txt", "T48hr_DON_H3K27me3_idr.05.annotation.txt", "T48hr_NoQ_H3K27me3_idr.05.annotation.txt")
T48hr_H3K27ac_peaks <- "T48hr_H3K27ac_annotation.txt"
T48hr_H3K27me3_peaks <- "T48hr_H3K27me3_annotation.txt"


T48hr_ac_pro_peaks <- lapply(T48hr_ac_peaks, function(x){get_ProPeak_genes(x)})
T48hr_me3_pro_peaks <- lapply(T48hr_me3_peaks, function(x){get_ProPeak_genes(x)})


T48hr_H3K27ac <- read.delim2(file = T48hr_H3K27ac_peaks)
T48hr_H3K27ac_pro <- T48hr_H3K27ac[abs(T48hr_H3K27ac$Distance.to.TSS) < 2000, ]
T48hr_H3K27me3 <- read.delim2(file = T48hr_H3K27me3_peaks)
T48hr_H3K27me3_pro <- T48hr_H3K27me3[abs(T48hr_H3K27me3$Distance.to.TSS) < 2000, ]


cpm <- read.csv(file = "CD8_RNAseq_cpm.csv")

if(FALSE){
  cpm_T48hr_CTL_ac <- cpm[cpm$symbol %in% T48hr_ac_pro_peaks[[1]]$Gene.Name, c(1:7, 44:47)]
  cpm_T48hr_CTL_me3 <- cpm[cpm$symbol %in% T48hr_me3_pro_peaks[[1]]$Gene.Name, c(1:7, 44:47)]
  cpm_T48hr_CTL <- rbind(cpm_T48hr_CTL_ac, cpm_T48hr_CTL_me3)
  cpm_T48hr_CTL$epi <- c(rep("H3K27ac", length(cpm_T48hr_CTL_ac$entrez)), rep("H3K27me3", length(cpm_T48hr_CTL_me3$entrez)))
  cpm_T48hr_CTL %>% pivot_longer(cols = 8:11, names_to = "SAMPLE", values_to = "CPM") %>% ggplot(aes(x = SAMPLE, y = log10(CPM + 0.1), color = epi)) + geom_boxplot() + coord_flip()
  ggsave(filename = "CTL_ac_me3_genes.tiff", dpi = 300, width = 4, height = 4)

  cpm_T48hr_CB839_ac <- cpm[cpm$symbol %in% T48hr_ac_pro_peaks[[2]]$Gene.Name, c(1:7, 48:51)]
  cpm_T48hr_CB839_me3 <- cpm[cpm$symbol %in% T48hr_me3_pro_peaks[[2]]$Gene.Name, c(1:7, 48:51)]
  cpm_T48hr_CB839 <- rbind(cpm_T48hr_CB839_ac, cpm_T48hr_CB839_me3)
  cpm_T48hr_CB839$epi <- c(rep("H3K27ac", length(cpm_T48hr_CB839_ac$entrez)), rep("H3K27me3", length(cpm_T48hr_CB839_me3$entrez)))
  cpm_T48hr_CB839 %>% pivot_longer(cols = 8:11, names_to = "SAMPLE", values_to = "CPM") %>% ggplot(aes(x = SAMPLE, y = log10(CPM + 0.1), color = epi)) + geom_boxplot() + coord_flip()
  ggsave(filename = "CB839_ac_me3_genes.tiff", dpi = 300, width = 4, height = 4)

  cpm_T48hr_DON_ac <- cpm[cpm$symbol %in% T48hr_ac_pro_peaks[[3]]$Gene.Name, c(1:7, 52:55)]
  cpm_T48hr_DON_me3 <- cpm[cpm$symbol %in% T48hr_me3_pro_peaks[[3]]$Gene.Name, c(1:7, 52:55)]
  cpm_T48hr_DON <- rbind(cpm_T48hr_DON_ac, cpm_T48hr_DON_me3)
  cpm_T48hr_DON$epi <- c(rep("H3K27ac", length(cpm_T48hr_DON_ac$entrez)), rep("H3K27me3", length(cpm_T48hr_DON_me3$entrez)))
  cpm_T48hr_DON %>% pivot_longer(cols = 8:11, names_to = "SAMPLE", values_to = "CPM") %>% ggplot(aes(x = SAMPLE, y = log10(CPM + 0.1), color = epi)) + geom_boxplot() + coord_flip()
  ggsave(filename = "DON_ac_me3_genes.tiff", dpi = 300, width = 4, height = 4)


  cpm_T48hr_NoQ_ac <- cpm[cpm$symbol %in% T48hr_ac_pro_peaks[[4]]$Gene.Name, c(1:7, 56:59)]
  cpm_T48hr_NoQ_me3 <- cpm[cpm$symbol %in% T48hr_me3_pro_peaks[[4]]$Gene.Name, c(1:7, 56:59)]
  cpm_T48hr_NoQ <- rbind(cpm_T48hr_NoQ_ac, cpm_T48hr_NoQ_me3)
  cpm_T48hr_NoQ$epi <- c(rep("H3K27ac", length(cpm_T48hr_NoQ_ac$entrez)), rep("H3K27me3", length(cpm_T48hr_NoQ_me3$entrez)))
  cpm_T48hr_NoQ %>% pivot_longer(cols = 8:11, names_to = "SAMPLE", values_to = "CPM") %>% ggplot(aes(x = SAMPLE, y = log10(CPM + 0.1), color = epi)) + geom_boxplot() + coord_flip()
  ggsave(filename = "NoQ_ac_me3_genes.tiff", dpi = 300, width = 4, height = 4)
}

cpm_mean <- cpm %>% mutate(T48hr_CTL = (X7342.MM.37 + X7342.MM.38 + X7342.MM.39 + X7342.MM.40)/4, T48hr_CB839 = c(X7342.MM.41 + X7342.MM.42 + X7342.MM.43 + X7342.MM.44)/4, T48hr_DON = (X7342.MM.45 + X7342.MM.46 + X7342.MM.47 + X7342.MM.48)/4, T48hr_NoQ = c(X7342.MM.49 + X7342.MM.50 + X7342.MM.51 + X7342.MM.52)/4)

cpm_mean_ac <- cpm_mean[cpm_mean$symbol %in% T48hr_H3K27ac_pro$Gene.Name, c(1:7, 60:63)]
wilcox.test(cpm_mean_ac$T48hr_CTL, cpm_mean_ac$T48hr_DON, paired = TRUE)
# cpm_mean_ac %>% write.csv(file = "H3K27ac_genes.csv")
cpm_mean_ac %>% pivot_longer(cols = 8:11, names_to = "TREATMENT", values_to = "CPM") %>% ggplot(aes(x = TREATMENT, y = log10(CPM + 0.1))) + geom_boxplot() + ggtitle("H3K27ac genes expression level")
ggsave(filename = "H3K27ac_genes_expression.tiff", dpi = 300, width = 4, height = 4)


cpm_mean_me3 <- cpm_mean[cpm_mean$symbol %in% T48hr_H3K27me3_pro$Gene.Name, c(1:7, 60:63)]
wilcox.test(cpm_mean_me3$T48hr_CTL, cpm_mean_me3$T48hr_DON, paired = TRUE)
# cpm_mean_me3 %>% write.csv(file = "H3K27me_genes.csv")
cpm_mean_me3 %>% pivot_longer(cols = 8:11, names_to = "TREATMENT", values_to = "CPM") %>% ggplot(aes(x = TREATMENT, y = log10(CPM + 0.1))) +geom_boxplot() + ggtitle("H3K27me3 genes expression level")
ggsave(filename = "H3K27me3_genes_expression.tiff", dpi = 300, width = 4, height = 4)


# GSVA analysis ---------------------------------
library(GSVA)
library(msigdbr)
library(pheatmap)

C7_immu_gene_sets = msigdbr(species = "mouse", category = "C7", subcategory = "IMMUNESIGDB")
C7_immu_gs = split(x = C7_immu_gene_sets$entrez_gene, f = C7_immu_gene_sets$gs_name)

C8_gene_sets = msigdbr(species = "mouse", category = "C8")
C8_gs = split(x = C8_gene_sets$entrez_gene, f = C8_gene_sets$gs_name)

C2_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG")
C2_gs = split(x = C2_gene_sets$entrez_gene, f = C2_gene_sets$gs_name)

C2_REACTOME_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:REACTOME")
C2_REACTOME_gs = split(x = C2_REACTOME_gene_sets$entrez_gene, f = C2_REACTOME_gene_sets$gs_name)

C1_gene_sets = msigdbr(species = "mouse", category = "C1")
C1_gs = split(x = C1_gene_sets$entrez_gene, f = C1_gene_sets$gs_name)

logNormCount <- cpm(y, log = TRUE)


col_ann <- data.frame(genotype = factor(metaData$group, levels = c("Vehicle.0hr",
                                                                   "Vehicle.4hr", "CB839.4hr", "DON.4hr", "No_Q.4hr",
                                                                   "Vehicle.24hr", "CB839.24hr", "DON.24hr", "No_Q.24hr",
                                                                   "Vehicle.48hr", "CB839.48hr", "DON.48hr", "No_Q.48hr")))
rownames(col_ann) <- colnames(logNormCount)

C8_en <- gsva(logNormCount, gset.idx.list = C8_gs)
KEGG_en <- gsva(logNormCount, gset.idx.list = C2_gs)
REACTOME_en <- gsva(logNormCount, gset.idx.list = C2_REACTOME_gs)
C1_en <- gsva(logNormCount, gset.idx.list = C1_gs)

#CD8_subset <- gs_en[names(C7_immu_gs) %in% names(C7_immu_gs)[grep("CD8", names(C7_immu_gs))],]

#pheatmap(CD8_subset, annotation_col = col_ann, cluster_cols = FALSE, show_rownames = TRUE)
pheatmap(C8_en, annotation_col = col_ann, cluster_cols = FALSE, show_rownames = TRUE)

pheatmap(KEGG_en, annotation_col = col_ann, cluster_cols = FALSE, show_rownames = TRUE)
pheatmap(KEGG_en, annotation_col = col_ann, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE)

pheatmap(C1_en, annotation_col = col_ann, cluster_cols = FALSE, show_rownames = TRUE)

pheatmap(REACTOME_en, annotation_col = col_ann, cluster_cols = FALSE, show_rownames = TRUE)

# for TCseq analysis ---------------------------------------------------------------
require(TCseq)
require(rtracklayer)
data("CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES")
data("CD8_RNAseq_GLUTAMINE_METABOLISM_COUNTS")
data("CD8_RNAseq_GLUTAMINE_METABOLISM_FEATURE_COUNT")


genes_interval <- import.gff2(con = "gencode.vM25.primary_assembly.annotation.gtf")
genes_interval <- genes_interval[mcols(genes_interval)$type %in% "gene"]
genes_interval <- data.frame(genes_interval)
genes_feature <- genes_interval[, c(1,2,3,5,10)]
colnames(genes_feature) <- c("chr", "start", "end", "strand", "id")

control_idx <- c(1:4, 5:8, 21:24, 37:40)

ctx <- as.matrix(CD8_RNAseq_GLUTAMINE_METABOLISM_FEATURE_COUNT[, 7:58])
rownames(ctx) <- CD8_RNAseq_GLUTAMINE_METABOLISM_FEATURE_COUNT$entrez

tca <- TCA(design = CD8_RNAseq_GLUTAMINE_METABOLISM_SAMPLES[control_idx,], genomicFeature = genes_feature, counts = ctx[, control_idx])
#tca <- countReads(tca, dir = "/home/xiang/DNA_seq/CD8_RNAseq/aligned")
tca <- DBanalysis(tca, filter.type = "raw", filter.value = 10, samplePassfilter = 4)
tca <- timecourseTable(tca, value = "expression", norm.method = "rpkm", filter = TRUE)

tca <- timeclust(tca, algo = "cm", k = 6)


#counts(tca)
