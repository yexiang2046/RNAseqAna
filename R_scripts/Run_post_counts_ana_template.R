require(devtools)
devtools::install_version("dbplyr", version = "2.3.4")

require(dplyr) 
require(readr)
require(edgeR)
require(clusterProfiler)
require(fgsea)
require(msigdbr)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(Organism.dplyr)
require(ggfortify) # for PCA plot
require(cowplot)
require(ggrepel)
require(tidyplots)
require(tibble)
require(factoextra)
require(VennDiagram)


load_all()


metaData <- readxl::read_excel(path = "SAMPLE_SHEET/12630-RS_project_summary.xlsx")
metaData$group <- paste0(metaData$Celltype, "_", metaData$genotype, "_", metaData$Temp)

metaData$SampleId <- metaData$Filename
metaData <- data.frame(metaData)

metaData %>% write.table("metaData_12630-RS.txt", sep = "\t", row.names = FALSE)

counts <- read.delim(file = "counts.txt", header = TRUE, skip = 1)

step1_get_counts_output <- step1_get_counts(featurecount_file = "counts.txt",
                                            type = "command_line",
                                            col_names = c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                                                        "12630.RS.019_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.023_S1_L005_RAligned.sortedByCoord.out.bam", 
                                                        "12630.RS.029_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.010_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.017_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.011_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.020_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.028_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.024_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.001_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.009_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.002_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.006_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.004_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.027_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.033_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.032_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.012_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.015_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.031_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.022_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.039_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.025_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.030_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.013_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.005_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.026_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.003_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.038_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.008_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.014_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.035_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.021_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.034_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.036_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.037_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.040_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.016_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.007_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.018_S1_L005_RAligned.sortedByCoord.out.bam"),
                                            counts_col = c("12630.RS.019_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.023_S1_L005_RAligned.sortedByCoord.out.bam", 
                                                        "12630.RS.029_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.010_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.017_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.011_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.020_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.028_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.024_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.001_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.009_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.002_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.006_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.004_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.027_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.033_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.032_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.012_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.015_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.031_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.022_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.039_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.025_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.030_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.013_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.005_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.026_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.003_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.038_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.008_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.014_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.035_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.021_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.034_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.036_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.037_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.040_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.016_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.007_S1_L005_RAligned.sortedByCoord.out.bam",
                                                        "12630.RS.018_S1_L005_RAligned.sortedByCoord.out.bam"))

Group <- factor(metaData$Group)
design <- model.matrix(~0 + Group)
colnames(design) <- make.names(levels(Group))
pairwise_contrasts <- combn(levels(Group), 2, function(x) {
  makeContrasts(paste0(make.names(x[1]), " - ", make.names(x[2])), levels = design)
}, simplify = FALSE)

con1 <- makeContrasts(GroupiTreg_cGASko_37 - GroupiTreg_cGASko_39, levels = design)
con2 <- makeContrasts(groupNSUN2_KO - groupCTR, levels = design)
con3 <- makeContrasts(groupNSUN2_KO - groupCMV_NSUN2, levels = design)

# need to be generalized
con <- list(con1, con2, con3)
step2_edgeR_analysis_output <- step2_edgeR_analysis(counts = step1_get_counts_output, meta = metaData,
                                                    sample_column = "SampleId",
                                                    groups = "group",
                                                    edgeR_file = "NSUN2_edgeR.rds",
                                                    cpm_file = "NSUN2_edgeR_cpm.csv",
                                                    contrast_list = con,
                                                    logFC = 1.0, padj = 0.05)

step3_prep_annotation_output <- step3_prep_annotation(sp = "human",
                                                      IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
                                                      info_types = c("entrez", "symbol", "ensembl"),
                                                      IdType = "ensembl")


step4_annotate_output_LANAsg_C1 <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
                                               annotation = step3_prep_annotation_output,
                                               Id_col = "ensembl",
                                               file_name = "BCBL1_LANAsg_9_12")

step4_annotate_output_LANAsg_C2 <- step4_annotate(qlf = step2_edgeR_analysis_output[[2]],
                                                  annotation = step3_prep_annotation_output,
                                                  Id_col = "ensembl",
                                                  file_name = "BCBL1_LANAsg_10_12")

keggC1 <- step5_export_KEGG_terms(edgeR_tb = step4_annotate_output_LANAsg_C1)
dotplot(keggC1)

keggC2 <- step5_export_KEGG_terms(edgeR_tb = step4_annotate_output_LANAsg_C2)
dotplot(keggC2)


# replicate consistency

gene_sets_Up <- list(repA = step4_annotate_output_LANAsg_C1[step4_annotate_output_LANAsg_C1$threshold == "Up", ]$ensembl,
           repB = step4_annotate_output_LANAsg_C2[step4_annotate_output_LANAsg_C2$threshold == "Up", ]$ensembl)

venn.diagram(x = gene_sets_Up, filename = "BCBL1_LANA_KD_replicateUpGenes_vennD.tiff")

gene_sets_Down <- list(repA = step4_annotate_output_LANAsg_C1[step4_annotate_output_LANAsg_C1$threshold == "Down", ]$ensembl,
                     repB = step4_annotate_output_LANAsg_C2[step4_annotate_output_LANAsg_C2$threshold == "Down", ]$ensembl)

venn.diagram(x = gene_sets_Down, filename = "BCBL1_LANA_KD_replicateDownGenes_vennD.tiff")


# combine two guide combinations --------------------
group <- factor(metaData[, "Treatment_agg"])
design <- model.matrix(~0+group)
con1 <- makeContrasts(groupLANAsg_KD - groupNT, levels = design)


# need to be generalized
con <- list(con1)
step2_edgeR_analysis_output <- step2_edgeR_analysis(counts = step1_get_counts_output, meta = metaData,
                                                    sample_column = "sampleID",
                                                    groups = "Treatment_agg",
                                                    edgeR_file = "BCBL1_LANA_KD_edgeR.rds",
                                                    cpm_file = "BCBL1_LANA_KD_edgeR_cpm.csv",
                                                    contrast_list = con,
                                                    logFC = 1.0, padj = 0.01)

step3_prep_annotation_output <- step3_prep_annotation(sp = "human",
                                                      IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
                                                      info_types = c("entrez", "symbol", "ensembl"),
                                                      IdType = "ensembl")


step4_annotate_output_LANAsg_KD <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
                                                  annotation = step3_prep_annotation_output,
                                                  Id_col = "ensembl",
                                                  file_name = "BCBL1_LANAsg_KD_4reps")

sig_degs <- step4_annotate_output_LANAsg_KD[step4_annotate_output_LANAsg_KD$padj < 0.05 & abs(step4_annotate_output_LANAsg_KD$logFC) > 1,]

sig_label <- step4_annotate_output_LANAsg_KD[step4_annotate_output_LANAsg_KD$padj < 1e-5 & abs(step4_annotate_output_LANAsg_KD$logFC) >5,]

step4_annotate_output_LANAsg_KD %>%
  mutate(significance = factor(c(ifelse((abs(logFC) > 1 & padj < 0.05), "Yes", "No")), levels = c("Yes", "No"))) %>%
  ggplot(aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(color = significance), size = 0.2) +
  scale_color_manual(name = "significant", values = c("red", "black")) +
  theme_classic() #+
#  geom_label_repel(data = step4_annotate_output_LANAsg_KD %>%
#                     filter(abs(logFC) > 5 & -log10(padj) > 5), aes(label = symbol), max.overlaps = 20)



keggLANA_KD <- step5_export_KEGG_terms(edgeR_tb = step4_annotate_output_LANAsg_KD)
dotplot(keggLANA_KD)



