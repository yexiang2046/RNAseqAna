
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

  # samples for Renca cells
  idx_use <- c(1:12)

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


  metaData <- metaData[idx_use,]

  group <- factor(paste(metaData$Genotype, metaData$Tissue, sep = "_"),
                  levels = unique(paste(metaData$Genotype, metaData$Tissue, sep = "_")))

  metaData <- cbind(metaData, group)
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  #con1 <- makeContrasts(Vhl_KO_Renca - Vhl_WT_Renca, levels = design)
  #con2 <- makeContrasts(Vhl_KO_Tumor - Vhl_WT_Tumor, levels = design)
  #con3 <- makeContrasts(Vhl_KO_MDSC - Vhl_WT_MDSC, levels = design)
  #con4 <- makeContrasts(Vhl_KO_PMN - Vhl_WT_PMN, levels = design)
  #con5 <- makeContrasts(Vhl_KO_F480lo_CD206lo - Vhl_WT_F480lo_CD206lo, levels = design)
  #con6 <- makeContrasts(Vhl_KO_F480hi_CD206hi - Vhl_WT_F480hi_CD206hi, levels = design)

  con7 <- makeContrasts(Vhl_KO1_Renca - Vhl_WT1_Renca, levels = design)
  con8 <- makeContrasts(Vhl_KO1_Renca - Vhl_WT3_Renca, levels = design)
  con9 <- makeContrasts(Vhl_KO2_Renca - Vhl_WT1_Renca, levels = design)
  con10 <- makeContrasts(Vhl_KO2_Renca - Vhl_WT3_Renca, levels = design)
  con11 <- makeContrasts(Vhl_KO2_Renca - Vhl_KO1_Renca, levels = design)
  con12 <- makeContrasts(Vhl_WT3_Renca - Vhl_WT1_Renca, levels = design)

  # need to be generalized
  # con <- list(con2, con3, con4, con5, con6, con7, con8, con9, con10, con11, con12)
  con <- list(con7, con8, con9, con10, con11, con12)



  counts_used <- step1_get_counts_output[, idx_use]
  counts_used <- counts_used[, c(4:12, 1:3)]

  # reorder counts columns
  # counts_used <- counts_used[, c(4:12, 1:3, 13:50)]
  step2_edgeR_analysis_output <- step2_edgeR_analysis(counts = counts_used, meta = metaData_used,
                                                      sample_column = "SampleID",
                                                      groups = "group",
                                                      edgeR_file = "Renca_cells_edgeR.rds",
                                                      cpm_file = "Renca_cells_edgeR_cpm.csv",
                                                      contrast_list = con,
                                                      logFC = 1.0, padj = 0.01)

  step3_prep_annotation_output <- step3_prep_annotation(sp = "mouse",
                                                        IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
                                                        info_types = c("entrez", "symbol", "ensembl"),
                                                        IdType = "ensembl")




  step4_annotate_output_Vhl_KO1_Renca_vs_Vhl_WT1_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
                                                                         annotation = step3_prep_annotation_output,
                                                                         Id_col = "ensembl",
                                                                         file_name = "Vhl_KO1_Renca_vs_Vhl_WT1_Renca_DEGs")
  step4_annotate_output_Vhl_KO1_Renca_vs_Vhl_WT2_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[2]],
                                                                         annotation = step3_prep_annotation_output,
                                                                         Id_col = "ensembl",
                                                                         file_name = "Vhl_KO1_Renca_vs_Vhl_WT2_Renca_DEGs")
  step4_annotate_output_Vhl_KO2_Renca_vs_Vhl_WT1_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[3]],
                                                                         annotation = step3_prep_annotation_output,
                                                                         Id_col = "ensembl",
                                                                         file_name = "Vhl_KO2_Renca_vs_Vhl_WT1_Renca_DEGs")
  step4_annotate_output_Vhl_KO2_Renca_vs_Vhl_WT2_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[4]],
                                                                         annotation = step3_prep_annotation_output,
                                                                         Id_col = "ensembl",
                                                                         file_name = "Vhl_KO2_Renca_vs_Vhl_WT2_Renca_DEGs")
  step4_annotate_output_Vhl_KO2_Renca_vs_Vhl_KO1_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[5]],
                                                                         annotation = step3_prep_annotation_output,
                                                                         Id_col = "ensembl",
                                                                         file_name = "Vhl_KO2_Renca_vs_Vhl_KO1_Renca_DEGs")
  step4_annotate_output_Vhl_WT2_Renca_vs_Vhl_WT1_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[6]],
                                                                         annotation = step3_prep_annotation_output,
                                                                         Id_col = "ensembl",
                                                                         file_name = "Vhl_WT2_Renca_vs_Vhl_WT1_Renca_DEGs")
}


y <- readRDS("Renca_cells_edgeR.rds")




colors <- rep(c("gray", "red", "black", "darkred"), 3)

points <- rep(1:4, 3)
plotMDS(y, pch = points, col = colors, gene.selection = "pairwise")
legend("topleft", legend=levels(metaData$group),
       pch=unique(points), col=unique(colors), ncol=3, cex = 0.8)




# export cpm with gene annotations ---------------------------------
library(Organism.dplyr)
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
cts %>% write_csv(file = "Renca_cells_cpm.csv")


# gene expression heatmap --------------
library(pheatmap)
mx <- as.matrix(cts[, 2:13])
rownames(mx) <- rownames(cts)
colnames(mx) <- colnames(cts)[2:13]



col_ann <- data.frame(genotype = factor(metaData$Genotype, levels = c("Vhl_WT1", "Vhl_WT3", "Vhl_KO1", "Vhl_KO2")))
rownames(col_ann) <- colnames(cts)[2:13]

pheatmap(mx, scale = "row",  annotation_col = col_ann, show_rownames = FALSE)


# GSVA analysis on Renca cells --------------------------------
library(GSVA)
library(msigdbr)
library(pheatmap)



H_gene_sets = msigdbr(species = "mouse", category = "H")
H_gs = split(x = H_gene_sets$ensembl_gene, f = H_gene_sets$gs_name)

Renca_cells_cpm <- as.matrix(cts[, 2:13])
rownames(Renca_cells_cpm) <- cts$ensembl

col_ann <- data.frame(genotype = factor(metaData$Genotype[1:12],levels = c("Vhl_WT1", "Vhl_KO1", "Vhl_WT3", "Vhl_KO2")))
rownames(col_ann) <- colnames(Renca_cells_cpm)

gs_en <- gsva(Renca_cells_cpm, gset.idx.list = H_gs)
gs_en <- gs_en[, c(1,5,9,3,7,11,2,6,10,4,8,12)]
pheatmap(gs_en, annotation_col = col_ann, cluster_cols = FALSE)


