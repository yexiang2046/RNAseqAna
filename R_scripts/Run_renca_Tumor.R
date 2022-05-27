
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



group <- factor(paste(metaData$Genotype, metaData$Tissue, sep = "_"),
                levels = unique(paste(metaData$Genotype, metaData$Tissue, sep = "_")))

metaData <- cbind(metaData, group)
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

# filter out failed samples
idx_use <- c(21, 23:29, 31:52)
metaData_used <- metaData[idx_use,]
counts_used <- step1_get_counts_output[, idx_use]

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



step4_annotate_output_Tumor <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
                                              annotation = step3_prep_annotation_output,
                                              Id_col = "ensembl",
                                              file_name = "Tumor_WT_KO_DEGs")
step4_annotate_output_MDSC <- step4_annotate(qlf = step2_edgeR_analysis_output[[2]],
                                             annotation = step3_prep_annotation_output,
                                             Id_col = "ensembl",
                                             file_name = "MDSC_WT_KO_DEGs")
step4_annotate_output_PMN <- step4_annotate(qlf = step2_edgeR_analysis_output[[3]],
                                            annotation = step3_prep_annotation_output,
                                            Id_col = "ensembl",
                                            file_name = "PMN_WT_KO_DEGs")
step4_annotate_output_F480lo_CD206lo <- step4_annotate(qlf = step2_edgeR_analysis_output[[4]],
                                                       annotation = step3_prep_annotation_output,
                                                       Id_col = "ensembl",
                                                       file_name = "F480lo_CD206lo_WT_KO_DEGs")
step4_annotate_output_F480hi_CD206hi <- step4_annotate(qlf = step2_edgeR_analysis_output[[5]],
                                                       annotation = step3_prep_annotation_output,
                                                       Id_col = "ensembl",
                                                       file_name = "F480hi_CD206hi_WT_KO_DEGs")

step4_annotate_output_Vhl_KO1_Renca_vs_Vhl_WT1_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[6]],
                                                                       annotation = step3_prep_annotation_output,
                                                                       Id_col = "ensembl",
                                                                       file_name = "Vhl_KO1_Renca_vs_Vhl_WT1_Renca_DEGs")
step4_annotate_output_Vhl_KO1_Renca_vs_Vhl_WT2_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[7]],
                                                                       annotation = step3_prep_annotation_output,
                                                                       Id_col = "ensembl",
                                                                       file_name = "Vhl_KO1_Renca_vs_Vhl_WT2_Renca_DEGs")
step4_annotate_output_Vhl_KO2_Renca_vs_Vhl_WT1_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[8]],
                                                                       annotation = step3_prep_annotation_output,
                                                                       Id_col = "ensembl",
                                                                       file_name = "Vhl_KO2_Renca_vs_Vhl_WT1_Renca_DEGs")
step4_annotate_output_Vhl_KO2_Renca_vs_Vhl_WT2_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[9]],
                                                                       annotation = step3_prep_annotation_output,
                                                                       Id_col = "ensembl",
                                                                       file_name = "Vhl_KO2_Renca_vs_Vhl_WT2_Renca_DEGs")
step4_annotate_output_Vhl_KO2_Renca_vs_Vhl_KO1_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[10]],
                                                                       annotation = step3_prep_annotation_output,
                                                                       Id_col = "ensembl",
                                                                       file_name = "Vhl_KO2_Renca_vs_Vhl_KO1_Renca_DEGs")
step4_annotate_output_Vhl_WT2_Renca_vs_Vhl_WT1_Renca <- step4_annotate(qlf = step2_edgeR_analysis_output[[5]],
                                                                       annotation = step3_prep_annotation_output,
                                                                       Id_col = "ensembl",
                                                                       file_name = "Vhl_WT2_Renca_vs_Vhl_WT1_Renca_DEGs")


y <- readRDS("myeloid_edgeR.rds")

idx_use <- c(1:21, 23:29, 31:52)



colors <- c(rep(c("gray", "red", "dark", "darkred"), 3), rep(c("green", "blue"), 4),
            rep(c("coral", "cyan", "deeppink", "deepskyblue",
                  "coral4", "cyan4", "deeppink4", "deepskyblue4"), 4)[c(-2,-10)])

points <- c(rep(1:4, 3), rep(c(5, 6), 4),
            rep(c(5, 6, 7, 8,
                  9, 10, 11, 12), 4)[c(-2,-10)])
plotMDS(y, pch = points[idx_use], col = colors[idx_use], gene.selection = "pairwise")
legend("top", legend=levels(metaData_used$group),
       pch=unique(points), col=unique(colors), ncol=3, cex = 0.8)

# myeloid only analysis------------------------------------------------------------------

# MDS plot ------------------------------
idx_myeloid <- c(21:50)
plotMDS(y[,idx_myeloid], pch = points[idx_myeloid], col = colors[idx_myeloid], gene.selection = "pairwise",
        dim.plot = c(1, 2))
legend("topright", legend=c("Vhl_WT_MDSC",
                            "Vhl_WT_PMN",
                            "Vhl_WT_F480lo_CD206lo",
                            "Vhl_WT_F480hi_CD206hi",
                            "Vhl_KO_MDSC",
                            "Vhl_KO_PMN",
                            "Vhl_KO_F480lo_CD206lo",
                            "Vhl_KO_F480hi_CD206hi"),
       pch=points[43:50], col=colors[43:50], ncol=3, cex = 0.60)

plotMDS(y[,idx_myeloid], pch = points[idx_myeloid], col = colors[idx_myeloid], gene.selection = "pairwise",
        dim.plot = c(1, 3))
legend("topright", legend=c("Vhl_WT_MDSC",
                            "Vhl_WT_PMN",
                            "Vhl_WT_F480lo_CD206lo",
                            "Vhl_WT_F480hi_CD206hi",
                            "Vhl_KO_MDSC",
                            "Vhl_KO_PMN",
                            "Vhl_KO_F480lo_CD206lo",
                            "Vhl_KO_F480hi_CD206hi"),
       pch=points[43:50], col=colors[43:50], ncol=3, cex = 0.60)





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
# cts %>% write_csv(file = "myeloid_cpm.csv")


# gene expression heatmap --------------
mye_cts <- t(scale(t(as.matrix(cts[, idx_myeloid]))))
rownames(mye_cts) <- rownames(cts)
colnames(mye_cts) <- colnames(cts)[idx_myeloid]

col_ann <- data.frame(genotype = metaData_used$Genotype[idx_myeloid],
                      cell_population = metaData_used$Tissue[idx_myeloid])
rownames(col_ann) <- colnames(mye_cts)

pheatmap(mye_cts, annotation_col = col_ann)



# GSVA analysis on Renca cells --------------------------------
library(GSVA)
library(msigdbr)
library(pheatmap)

c2KEGG_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG")
msigdbr_list = split(x = c2KEGG_gene_sets$ensembl_gene, f = c2KEGG_gene_sets$gs_name)

H_gene_sets = msigdbr(species = "mouse", category = "H")
H_gs = split(x = H_gene_sets$ensembl_gene, f = H_gene_sets$gs_name)

Renca_cells_cpm <- as.matrix(cts[, 2:13])
rownames(Renca_cells_cpm) <- cts$ensembl

col_ann <- data.frame(genotype = factor(metaData$Genotype[1:12],levels = c("Vhl_WT1", "Vhl_WT3", "Vhl_KO1", "Vhl_KO2")))
rownames(col_ann) <- colnames(Renca_cells_cpm)

gs_en <- gsva(Renca_cells_cpm, gset.idx.list = msigdbr_list)
pheatmap(gs_en, annotation_col = col_ann)


gs_en_H <- gsva(Renca_cells_cpm, gset.idx.list = H_gs, method = "gsva")
pheatmap(gs_en_H, annotation_col = col_ann)


# GSVA analysis on Myeloid cells ---------------------------

Myeloid_cells_cpm <- as.matrix(cts[, c(22:51)])
rownames(Myeloid_cells_cpm) <- cts$ensembl

col_ann <- data.frame(genotype = factor(metaData$Genotype[idx_use[21:50]],
                                        levels = c("Vhl_WT", "Vhl_KO")),
                      cell_subpopulation = factor(metaData$Tissue[idx_use[21:50]]))
rownames(col_ann) <- colnames(Myeloid_cells_cpm)
gs_mye_en <- gsva(Myeloid_cells_cpm, gset.idx.list = H_gs)
pheatmap(gs_mye_en, annotation_col = col_ann)



