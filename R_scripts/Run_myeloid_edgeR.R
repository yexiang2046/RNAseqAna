
require(dplyr)
require(readr)
require(readxl)
require(edgeR)
require(clusterProfiler)
require(fgsea)
require(msigdbr)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Mm.eg.db)
require(org.Hs.eg.db)
require(Organism.dplyr)
require(ggfortify) # for PCA plot
require(cowplot)
require(ggrepel)
require(tibble)
require(factoextra)
require(VennDiagram)


require(devtools)
load_all()

if(FALSE){

  d1 <- read_xlsx(path = "7247-MW_7247_Wolf_Internal_Nucleic_Acid_Extraction_NGS_Submission_Form.xlsx", skip = 28)
  d2 <- read_xlsx(path = "7525-MW_7525_Wolf_Internal_Nucleic_Acid_Extraction_NGS_Submission_Form.xlsx", skip = 68)


  metaData <- data.frame(FileNameBase = c(d1$`Sample IDs`, d2$...2),
                         SampleID = paste(c(rep("Renca", 12), rep("Tumor", 8), rep(c("MDSC", "PMN", "F480lo_CD206lo", "F480hi_CD206hi"), 8)),
                                          c(rep(c("WT1", "WT3", "KO1", "KO2"), 3),
                                            rep(c("Vhl_WT", "Vhl_KO"), 4),
                                            rep(c("Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO"), each = 4)), sep = "_"),
                         Genotype = c(rep(c("Vhl_WT1", "Vhl_WT3", "Vhl_KO1","Vhl_KO2"), 3),
                                      rep(c("Vhl_WT", "Vhl_KO"), 4),
                                      rep(c("Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO"), each = 4)),
                         Tissue = c(rep("Renca", 12), rep("Tumor", 8), rep(c("MDSC", "PMN", "F480lo_CD206lo", "F480hi_CD206hi"), 8)))




  step1_get_counts_output <- step1_get_counts(featurecount_file = "myeloid_count_genes_vM25.txt",
                                              type = "command_line",
                                              col_names = c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                                                            metaData$FileNameBase[c(10:12, 1:9, 13:52)]),
                                              counts_col = metaData$FileNameBase[c(10:12, 1:9, 13:52)])


  # myeloid samples
  idx_use <- c(21, 23:29, 31:52)
  metaData_used <- metaData[idx_use,]
  counts_used <- step1_get_counts_output[, idx_use]


  group <- factor(paste(metaData_used$Genotype, metaData_used$Tissue, sep = "_"),
                  levels = unique(paste(metaData_used$Genotype, metaData_used$Tissue, sep = "_")))
  metaData_used <- cbind(metaData_used, group)

  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  #con1 <- makeContrasts(Vhl_KO_Renca - Vhl_WT_Renca, levels = design)
  #con2 <- makeContrasts(Vhl_KO_Tumor - Vhl_WT_Tumor, levels = design)
  con3 <- makeContrasts(Vhl_KO_MDSC - Vhl_WT_MDSC, levels = design)
  con4 <- makeContrasts(Vhl_KO_PMN - Vhl_WT_PMN, levels = design)
  con5 <- makeContrasts(Vhl_KO_F480lo_CD206lo - Vhl_WT_F480lo_CD206lo, levels = design)
  con6 <- makeContrasts(Vhl_KO_F480hi_CD206hi - Vhl_WT_F480hi_CD206hi, levels = design)

  #con7 <- makeContrasts(Vhl_KO1_Renca - Vhl_WT1_Renca, levels = design)
  #con8 <- makeContrasts(Vhl_KO1_Renca - Vhl_WT3_Renca, levels = design)
  #con9 <- makeContrasts(Vhl_KO2_Renca - Vhl_WT1_Renca, levels = design)
  #con10 <- makeContrasts(Vhl_KO2_Renca - Vhl_WT3_Renca, levels = design)
  #con11 <- makeContrasts(Vhl_KO2_Renca - Vhl_KO1_Renca, levels = design)
  #con12 <- makeContrasts(Vhl_WT3_Renca - Vhl_WT1_Renca, levels = design)

  # need to be generalized
  # con <- list(con2, con3, con4, con5, con6, con7, con8, con9, con10, con11, con12)
  con <- list(con3, con4, con5, con6)



  # reorder counts columns
  # counts_used <- counts_used[, c(4:12, 1:3, 13:50)]
  step2_edgeR_analysis_output <- step2_edgeR_analysis(counts = counts_used, meta = metaData_used,
                                                      sample_column = "SampleID",
                                                      groups = "group",
                                                      edgeR_file = "myeloid_edgeR.rds",
                                                      cpm_file = "myeloid_edgeR_cpm.csv",
                                                      contrast_list = con,
                                                      logFC = 1.0, padj = 0.01)

  step3_prep_annotation_output <- step3_prep_annotation(sp = "mouse",
                                                        IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
                                                        info_types = c("entrez", "symbol", "ensembl"),
                                                        IdType = "ensembl")




  step4_annotate_output_MDSC <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
                                               annotation = step3_prep_annotation_output,
                                               Id_col = "ensembl",
                                               file_name = "MDSC_WT_KO_DEGs")
  step4_annotate_output_PMN <- step4_annotate(qlf = step2_edgeR_analysis_output[[2]],
                                              annotation = step3_prep_annotation_output,
                                              Id_col = "ensembl",
                                              file_name = "PMN_WT_KO_DEGs")
  step4_annotate_output_F480lo_CD206lo <- step4_annotate(qlf = step2_edgeR_analysis_output[[3]],
                                                         annotation = step3_prep_annotation_output,
                                                         Id_col = "ensembl",
                                                         file_name = "F480lo_CD206lo_WT_KO_DEGs")
  step4_annotate_output_F480hi_CD206hi <- step4_annotate(qlf = step2_edgeR_analysis_output[[4]],
                                                         annotation = step3_prep_annotation_output,
                                                         Id_col = "ensembl",
                                                         file_name = "F480hi_CD206hi_WT_KO_DEGs")

}

require(dplyr)
require(readr)

require(clusterProfiler)
require(fgsea)
require(msigdbr)
require(org.Mm.eg.db)
require(org.Hs.eg.db)
require(Organism.dplyr)
require(ggfortify) # for PCA plot
require(cowplot)
require(ggrepel)
require(tibble)
# require(factoextra)
# require(VennDiagram)

require(readxl)
require(edgeR)
require(devtools)
load_all()

y <- readRDS("RDS/myeloid_edgeR.rds")

d1 <- read_xlsx(path = "SAMPLE_SHEET/7247-MW_7247_Wolf_Internal_Nucleic_Acid_Extraction_NGS_Submission_Form.xlsx", skip = 28)
d2 <- read_xlsx(path = "SAMPLE_SHEET/7525-MW_7525_Wolf_Internal_Nucleic_Acid_Extraction_NGS_Submission_Form.xlsx", skip = 68)


metaData <- data.frame(FileNameBase = c(d1$`Sample IDs`, d2$...2),
                       SampleID = paste(c(rep("Renca", 12), rep("Tumor", 8), rep(c("MDSC", "PMN", "F480lo_CD206lo", "F480hi_CD206hi"), 8)),
                                        c(rep(c("WT1", "WT3", "KO1", "KO2"), 3),
                                          rep(c("Vhl_WT", "Vhl_KO"), 4),
                                          rep(c("Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO"), each = 4)), sep = "_"),
                       Genotype = c(rep(c("Vhl_WT1", "Vhl_WT3", "Vhl_KO1","Vhl_KO2"), 3),
                                    rep(c("Vhl_WT", "Vhl_KO"), 4),
                                    rep(c("Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO", "Vhl_WT", "Vhl_KO"), each = 4)),
                       Tissue = c(rep("Renca", 12), rep("Tumor", 8), rep(c("MDSC", "PMN", "F480lo_CD206lo", "F480hi_CD206hi"), 8)))



idx_use <- c(21, 23:29, 31:52)

#use only TAM idx
idx_use_modified <- c(21,23,24,25,27,28,29,31,32,33,35,36,37,39,40,41,43,44,45,47,48,49,51,52)
metaData_used <- metaData[idx_use_modified,]
group <- factor(paste(metaData_used$Genotype, metaData_used$Tissue, sep = "_"),
                levels = unique(paste(metaData_used$Genotype, metaData_used$Tissue, sep = "_")))
metaData_used <- cbind(metaData_used, group)

colors <- c(rep(c("coral", "cyan", "deeppink", "deepskyblue",
                  "coral4", "cyan4", "deeppink4", "deepskyblue4"), 4)[c(-2,-10)])

points <- c(rep(c(5, 6, 7, 8,
                  9, 10, 11, 12), 4)[c(-2,-10)])
plotMDS(y, pch = points, col = colors, gene.selection = "pairwise")
legend("topright", legend=levels(metaData_used$group),
       pch=unique(points), col=unique(colors), ncol=3, cex = 0.6)






# heatmaps ---------------------------------
library(Organism.dplyr)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
cts <- cpm(y)
cts <- cts %>% data.frame() %>% rownames_to_column()
colnames(cts)[1] <- "ensembl"

src <- src_ucsc(organism = "mouse")
IdList <- sapply(cts$ensembl, function(x){strsplit(x, "\\.")}[[1]][[1]])

cts$ensembl <- IdList
sym <- select(src, keys = IdList, columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")

cts <- left_join(cts, sym, by = "ensembl")
cts <- cts[!duplicated(cts$ensembl),]
cts %>% write_csv(file = "myeloid_only_cpm.csv")



# gene expression heatmap --------------
library(pheatmap)
mx <- as.matrix(cts[, 2:31])
rownames(mx) <- rownames(cts)
colnames(mx) <- colnames(cts)[2:31]

col_order <- c(1,7,13,19,4,10,16,22,2,8,14,20,5,11,17,23,
3,9,15,21,6,12,18,24)

col_ann <- data.frame(genotype = factor(metaData_used$Genotype, levels = c("Vhl_WT", "Vhl_KO")),
                      cell_population = factor(metaData_used$Tissue, levels = c("MDSC", "PMN", "F480lo_CD206lo", "F480hi_CD206hi")))
rownames(col_ann) <- metaData[idx_use_modified, "FileNameBase"]

pheatmap(mx[, col_order], scale = "row",  annotation_col = col_ann, cluster_cols = FALSE, show_rownames = FALSE)








# GSVA analysis on Myeloid cells ---------------------------
library(GSVA)
library(msigdbr)
library(pheatmap)

h_gene_sets = msigdbr(species = "mouse", category = "H")

msigdbr_list = split(x = h_gene_sets$ensembl_gene, f = h_gene_sets$gs_name)

# Myeloid_cells_cpm <- as.matrix(cts[, c(2:31)])
Myeloid_cells_cpm <- as.matrix(cts[, c(2,3,4,5,6,7,9,10,11,13,14,15,17,18,19,21,22,23,25,26,27,29,30,31)])
rownames(Myeloid_cells_cpm) <- cts$ensembl

col_ann <- data.frame(genotype = factor(metaData$Genotype[idx_use_modified],
                                        levels = c("Vhl_WT", "Vhl_KO")),
                      cell_subpopulation = factor(metaData$Tissue[idx_use_modified],
                                                  levels = c("PMN", "MDSC", "F480lo_CD206lo", "F480hi_CD206hi")))
rownames(col_ann) <- colnames(Myeloid_cells_cpm)
gs_mye_en <- gsva(Myeloid_cells_cpm, gset.idx.list = msigdbr_list)

#col_order <- c(16, 24, 5, 12, 20,28,
#               1, 8,  15,  23,4, 11,19,27,
#               2,  9, 17, 25,  6, 13, 21, 29,
#               3, 10, 18,  26,  7, 14, 22, 30)
# pheatmap(gs_mye_en[, col_order], annotation_col = col_ann, cluster_cols = FALSE)

col_order <- c(1,7,13,19,4,10,16,22,2,8,14,20,5,11,17,23,3,9,15,21,6,12,18,24)
pheatmap(gs_mye_en[,col_order], annotation_col = col_ann[col_order,], cluster_cols = FALSE, show_colnames = FALSE)


# pearson correlation between samples
idx <- c(1:6,8:10,12:14,16:18,20:22,24:26,28:30)
colnames(cpm(y, log = TRUE))
rlog.norm.counts <- cpm(y, log = TRUE)
count_sample <- colnames(rlog.norm.counts)
names(count_sample) <- 1:30
rlog.norm.counts <- cpm(y, log = TRUE)[,idx]


col_ann <- data.frame(genotype = factor(metaData$Genotype[idx_use_modified],
                                        levels = c("Vhl_WT", "Vhl_KO")),
                      cell_subpopulation = factor(metaData$Tissue[idx_use_modified],
                                                  levels = c("PMN", "MDSC", "F480lo_CD206lo", "F480hi_CD206hi")))
col_ann <- data.frame(genotype = factor(metaData$Genotype[idx_use_modified],
                                        levels = c("Vhl_WT", "Vhl_KO")),
                      cell_subpopulation = factor(metaData$Tissue[idx_use_modified],
                                                  levels = c("PMN", "MDSC", "F480lo_CD206lo", "F480hi_CD206hi")))
rownames(col_ann) <- colnames(rlog.norm.counts)
col_order <- c(1,7,13,19,4,10,16,22,2,8,14,20,5,11,17,23,3,9,15,21,6,12,18,24)
corr_coeff <- cor(rlog.norm.counts[,col_order], method = "pearson") 
as.dist(1-corr_coeff, upper = TRUE) %>%
as.matrix %>%
pheatmap::pheatmap(., main = "Pearson correlation", cluster_cols = FALSE, annotation_col = col_ann)
