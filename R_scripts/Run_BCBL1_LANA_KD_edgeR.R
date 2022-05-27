
require(dplyr)
require(readr)
require(edgeR)
require(clusterProfiler)
require(fgsea)
require(msigdbr)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
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




metaData <- data.frame(sampleID = c(paste0(7505, sep = '_', 1:6)),
                       Treatment = c("NT", "LANAsg_C1", "LANAsg_C2", "NT", "LANAsg_C1", "LANAsg_C2"),
                       Treatment_agg = c("NT", "LANAsg_KD", "LANAsg_KD", "NT", "LANAsg_KD", "LANAsg_KD"))




step1_get_counts_output <- step1_get_counts(featurecount_file = "BCBL1_dSg_LANA_KD_RNAseq_counts.txt",
                                            type = "command_line",
                                            col_names = c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                                                          "NT_A", "LANAsg_C1_A", "LANAsg_C2_A",
                                                          "NT_B", "LANAsg_C1_B", "LANAsg_C2_B"),
                                            counts_col = c("NT_A", "LANAsg_C1_A", "LANAsg_C2_A",
                                                           "NT_B", "LANAsg_C1_B", "LANAsg_C2_B"))
## don't run this, combined replicates for DEGs ---------------------------------------
# group <- factor(metaData[, "Treatment"])
# design <- model.matrix(~0+group)
# con1 <- makeContrasts(groupLANAsg_C1 - groupNT, levels = design)
# con2 <- makeContrasts(groupLANAsg_C2 - groupNT, levels = design)
#
## need to be generalized
#con <- list(con1, con2)
#step2_edgeR_analysis_output <- step2_edgeR_analysis(counts = step1_get_counts_output, meta = metaData,
#                                                    sample_column = "sampleID",
#                                                    groups = "Treatment",
#                                                    edgeR_file = "BCBL1_LANA_KD_edgeR.rds",
#                                                    cpm_file = "BCBL1_LANA_KD_edgeR_cpm.csv",
#                                                    contrast_list = con,
#                                                    logFC = 1.0, padj = 0.01)
#
#step3_prep_annotation_output <- step3_prep_annotation(sp = "human",
#                                                      IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
#                                                      info_types = c("entrez", "symbol", "ensembl"),
#                                                      IdType = "ensembl")
#
#
#step4_annotate_output_LANAsg_C1 <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
#                                               annotation = step3_prep_annotation_output,
#                                               Id_col = "ensembl",
#                                               file_name = "BCBL1_LANAsg_9_12")
#
#step4_annotate_output_LANAsg_C2 <- step4_annotate(qlf = step2_edgeR_analysis_output[[2]],
#                                                  annotation = step3_prep_annotation_output,
#                                                  Id_col = "ensembl",
#                                                  file_name = "BCBL1_LANAsg_10_12")
#
#
#keggC1 <- step5_export_KEGG_terms(edgeR_tb = step4_annotate_output_LANAsg_C1)
#dotplot(keggC1)
#
#keggC2 <- step5_export_KEGG_terms(edgeR_tb = step4_annotate_output_LANAsg_C2)
#dotplot(keggC2)
#
#
## replicate consistency
#
#gene_sets_Up <- list(repA = step4_annotate_output_LANAsg_C1[step4_annotate_output_LANAsg_C1$threshold == "Up", ]$ensembl,
#           repB = step4_annotate_output_LANAsg_C2[step4_annotate_output_LANAsg_C2$threshold == "Up", ]$ensembl)
#
#venn.diagram(x = gene_sets_Up, filename = "BCBL1_LANA_KD_replicateUpGenes_vennD.tiff")
#
#gene_sets_Down <- list(repA = step4_annotate_output_LANAsg_C1[step4_annotate_output_LANAsg_C1$threshold == "Down", ]$ensembl,
#                     repB = step4_annotate_output_LANAsg_C2[step4_annotate_output_LANAsg_C2$threshold == "Down", ]$ensembl)
#
#venn.diagram(x = gene_sets_Down, filename = "BCBL1_LANA_KD_replicateDownGenes_vennD.tiff")


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



