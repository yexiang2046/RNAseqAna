
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


require(devtools)
load_all()




metaData <- data.frame(sampleID = c(paste0(8004, sep = '_', 1:16)),
                       Treatment = c("WT_ui_input", "WT_ui_input", "WT_in_input", "WT_in_input", "KO_ui_input",
                                     "KO_ui_input", "KO_in_input", "KO_in_input", "WT_ui_IP", "WT_ui_IP", "WT_in_IP",
                                     "WT_in_IP", "KO_ui_IP", "KO_ui_IP", "KO_in_IP", "KO_in_IP"),
                       SampleName = c("WT_ui_inputA", "WT_ui_inputB", "WT_in_inputA", "WT_in_inputB", "KO_ui_inputA",
                                      "KO_ui_inputB", "KO_in_inputA", "KO_in_inputB", "WT_ui_IPA", "WT_ui_IPB", "WT_in_IPA",
                                      "WT_in_IPB", "KO_ui_IPA", "KO_ui_IPB", "KO_in_IPA", "KO_in_IPB"))




step1_get_counts_output <- step1_get_counts(featurecount_file = "counts/featureCounts_GRCh38_R40_ctx.txt",
                                            type = "command_line",
                                            col_names = c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                                                          "WT_ui_inputA", "WT_ui_inputB", "WT_in_inputA", "WT_in_inputB", "KO_ui_inputA",
                                                          "KO_ui_inputB", "KO_in_inputA", "KO_in_inputB", "WT_ui_IPA", "WT_ui_IPB", "WT_in_IPA",
                                                          "WT_in_IPB", "KO_ui_IPA", "KO_ui_IPB", "KO_in_IPA", "KO_in_IPB"),
                                            counts_col = c("WT_ui_inputA", "WT_ui_inputB", "WT_in_inputA", "WT_in_inputB", "KO_ui_inputA",
                                                           "KO_ui_inputB", "KO_in_inputA", "KO_in_inputB", "WT_ui_IPA", "WT_ui_IPB", "WT_in_IPA",
                                                           "WT_in_IPB", "KO_ui_IPA", "KO_ui_IPB", "KO_in_IPA", "KO_in_IPB"))

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
group <- factor(metaData[, "Treatment"])
design <- model.matrix(~0+group)
con1 <- makeContrasts(groupWT_in_input - groupWT_ui_input, levels = design)
con2 <- makeContrasts(groupKO_in_input - groupKO_ui_input, levels = design)
con3 <- makeContrasts(groupKO_ui_input - groupWT_ui_input, levels = design)
con4 <- makeContrasts(groupKO_in_input - groupWT_in_input, levels = design)

con5 <- makeContrasts(groupWT_ui_IP - groupWT_ui_input, levels = design)
con6 <- makeContrasts(groupWT_in_IP - groupWT_in_input, levels = design)
con7 <- makeContrasts(groupKO_ui_IP - groupKO_ui_input, levels = design)
con8 <- makeContrasts(groupKO_in_IP - groupKO_in_input, levels = design)

con9 <- makeContrasts(groupKO_ui_IP - groupWT_ui_IP, levels = design)
con10 <- makeContrasts(groupKO_in_IP - groupWT_in_IP, levels = design)
con11 <- makeContrasts(groupWT_in_IP - groupWT_ui_IP, levels = design)
con12 <- makeContrasts(groupKO_in_IP - groupKO_ui_IP, levels = design)

# need to be generalized
con <- list(con1, con2, con3, con4, con5, con6, con7, con8, con9, con10, con11, con12)
step2_edgeR_analysis_output <- step2_edgeR_analysis(counts = step1_get_counts_output, meta = metaData,
                                                    sample_column = "sampleID",
                                                    groups = "Treatment",
                                                    edgeR_file = "A549_VACV_edgeR.rds",
                                                    cpm_file = "A549_VACV_cpm.csv",
                                                    contrast_list = con,
                                                    logFC = 1.0, padj = 0.01)

step3_prep_annotation_output <- step3_prep_annotation(sp = "human",
                                                      IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
                                                      info_types = c("entrez", "symbol", "ensembl"),
                                                      IdType = "ensembl")


step4_annotate_output_input_WT_in_ui <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
                                                  annotation = step3_prep_annotation_output,
                                                  Id_col = "ensembl",
                                                  file_name = "A549_VACV_input_WT_infect_vs_uninfect")
step4_annotate_output_input_KO_in_ui <- step4_annotate(qlf = step2_edgeR_analysis_output[[2]],
                                                  annotation = step3_prep_annotation_output,
                                                  Id_col = "ensembl",
                                                  file_name = "A549_VACV_input_KO_infect_vs_uninfect")
step4_annotate_output_input_ui_KO_WT <- step4_annotate(qlf = step2_edgeR_analysis_output[[3]],
                                                  annotation = step3_prep_annotation_output,
                                                  Id_col = "ensembl",
                                                  file_name = "A549_VACV_input_uninfect_KO_vs_WT")
step4_annotate_output_input_in_KO_WT <- step4_annotate(qlf = step2_edgeR_analysis_output[[4]],
                                                  annotation = step3_prep_annotation_output,
                                                  Id_col = "ensembl",
                                                  file_name = "A549_VACV_input_infect_KO_vs_WT")


step4_annotate_output_WT_uni_IP_vs_input <- step4_annotate(qlf = step2_edgeR_analysis_output[[5]],
                                                 annotation = step3_prep_annotation_output,
                                                 Id_col = "ensembl",
                                                 file_name = "A549_VACV_WT_uni_IP_vs_input")
step4_annotate_output_WT_inf_IP_vs_input <- step4_annotate(qlf = step2_edgeR_analysis_output[[6]],
                                                 annotation = step3_prep_annotation_output,
                                                 Id_col = "ensembl",
                                                 file_name = "A549_VACV_WT_inf_IP_vs_input")
step4_annotate_output_KO_uni_IP_vs_input <- step4_annotate(qlf = step2_edgeR_analysis_output[[7]],
                                                 annotation = step3_prep_annotation_output,
                                                 Id_col = "ensembl",
                                                 file_name = "A549_VACV_KO_uni_IP_vs_input")
step4_annotate_output_KO_inf_IP_vs_input <- step4_annotate(qlf = step2_edgeR_analysis_output[[8]],
                                                 annotation = step3_prep_annotation_output,
                                                 Id_col = "ensembl",
                                                 file_name = "A549_VACV_KO_inf_IP_vs_input")


step4_annotate_output_IP_uni_KO_vs_WT <- step4_annotate(qlf = step2_edgeR_analysis_output[[9]],
                                                       annotation = step3_prep_annotation_output,
                                                       Id_col = "ensembl",
                                                       file_name = "A549_VACV_IP_uni_KO_vs_WT")
step4_annotate_output_IP_inf_KO_WT <- step4_annotate(qlf = step2_edgeR_analysis_output[[10]],
                                                       annotation = step3_prep_annotation_output,
                                                       Id_col = "ensembl",
                                                       file_name = "A549_VACV_IP_inf_KO_vs_WT")
step4_annotate_output_IP_WT_inf_vs_uni <- step4_annotate(qlf = step2_edgeR_analysis_output[[11]],
                                                       annotation = step3_prep_annotation_output,
                                                       Id_col = "ensembl",
                                                       file_name = "A549_VACV_IP_WT_inf_vs_uni")
step4_annotate_output_IP_KO_inf_vs_uni <- step4_annotate(qlf = step2_edgeR_analysis_output[[12]],
                                                       annotation = step3_prep_annotation_output,
                                                       Id_col = "ensembl",
                                                       file_name = "A549_VACV_IP_KO_inf_vs_uni")




step4_annotate_output_input_ui_KO_WT %>%
  mutate(significance = factor(c(ifelse((abs(logFC) > 1 & padj < 0.05), "Yes", "No")), levels = c("Yes", "No"))) %>%
  ggplot(aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(color = significance), size = 0.2) +
  scale_color_manual(name = "significant", values = c("red", "black")) +
  theme_classic() #+
#  geom_label_repel(data = step4_annotate_output_LANAsg_KD %>%
#                     filter(abs(logFC) > 5 & -log10(padj) > 5), aes(label = symbol), max.overlaps = 20)



uni_WT_KO <- data.frame(uni_WT_IP_enrich_fold = step4_annotate_output_WT_uni_IP_vs_input$logFC,
                        uni_KO_IP_enrich_fold = step4_annotate_output_KO_uni_IP_vs_input$logFC)
uni_WT_KO %>% ggplot(aes(x = uni_WT_IP_enrich_fold, y = uni_KO_IP_enrich_fold)) + geom_point() +
  xlim(-15, 15)+ ylim(-15, 15)



inf_WT_KO <- data.frame(inf_WT_IP_enrich_fold = step4_annotate_output_WT_inf_IP_vs_input$logFC,
                        inf_KO_IP_enrich_fold = step4_annotate_output_KO_inf_IP_vs_input$logFC)
inf_WT_KO %>% ggplot(aes(x = inf_WT_IP_enrich_fold, y = inf_KO_IP_enrich_fold)) + geom_point()+
  xlim(-15, 15)+ ylim(-15, 15)


WT_inf_uni <- data.frame(WT_uni_IP_enrich_fold = step4_annotate_output_WT_uni_IP_vs_input$logFC,
                         WT_inf_IP_enrich_fold = step4_annotate_output_WT_inf_IP_vs_input$logFC)
WT_inf_uni %>% ggplot(aes(x = WT_uni_IP_enrich_fold, y = WT_inf_IP_enrich_fold)) + geom_point()+
  xlim(-15, 15)+ ylim(-15, 15)

