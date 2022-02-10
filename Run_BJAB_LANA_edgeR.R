
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
require(tibble)




metaData <- data.frame(sampleID = c(paste0(7183, sep = '_', 1:4)),
                       Treatment = c("GFP", "GFP", "LANA", "LANA"))




step1_get_counts_output <- step1_get_counts(featurecount_file = "BJABs_LANA_RNA_count.txt", type = "command_line")
group <- factor(metaData[, "Treatment"])
design <- model.matrix(~0+group)
con1 <- makeContrasts(groupLANA - groupGFP, levels = design) # need to be generalized
con <- list(con1)
step2_edgeR_analysis_output <- step2_edgeR_analysis(counts = step1_get_counts_output, meta = metaData,
                                                    sample_column = "sampleID",
                                                    groups = "Treatment",
                                                    contrast_list = con,
                                                    logFC = 1.0, padj = 0.01)

step3_prep_annotation_output <- step3_prep_annotation(sp = "human",
                                                      IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
                                                      info_types = c("entrez", "symbol", "ensembl"),
                                                      IdType = "ensembl")


step4_annotate_output_LANA_GFP <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
                                               annotation = step3_prep_annotation_output,
                                               Id_col = "ensembl",
                                               file_name = "BJABs_LANA_vs_GFP")




sig_degs <- step4_annotate_output_LANA_GFP[step4_annotate_output_LANA_GFP$padj < 0.05 & abs(step4_annotate_output_LANA_GFP$logFC) > 1,]

sig_label <- step4_annotate_output_LANA_GFP[step4_annotate_output_LANA_GFP$padj < 1e-5 & abs(step4_annotate_output_LANA_GFP$logFC) >5,]

step4_annotate_output_LANA_GFP %>%
  mutate(significance = factor(c(ifelse((abs(logFC) > 1 & padj < 0.05), "Yes", "No")), levels = c("Yes", "No"))) %>%
  ggplot(aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(color = significance), size = 0.2) +
  scale_color_manual(name = "significant", values = c("red", "black")) +
  theme_classic() +
  geom_label_repel(data = step4_annotate_output_LANA_GFP %>%
                     filter(abs(logFC) > 2 & -log10(padj) > 5), aes(label = symbol), max.overlaps = 50)

BJABsLANA_OE <- step5_export_KEGG_terms(edgeR_tb = step4_annotate_output_LANA_GFP)
dotplot(BJABsLANA_OE)

LANA_Up_genes <-
  step4_annotate_output_LANA_GFP[step4_annotate_output_LANA_GFP$padj < 0.05 &
                                   step4_annotate_output_LANA_GFP$logFC > 1,]

LANA_Down_genes <- step4_annotate_output_LANA_GFP[step4_annotate_output_LANA_GFP$padj < 0.05 &
                                                    step4_annotate_output_LANA_GFP$logFC < 1,]


BJABsLANA_OE_up <- step5_export_KEGG_terms(edgeR_tb = LANA_Up_genes)
dotplot(BJABsLANA_OE_up)

BJABsLANA_OE_down <- step5_export_KEGG_terms(edgeR_tb = LANA_Down_genes)
dotplot(BJABsLANA_OE_down)

kegg_gene_sets = msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
m2tg <- kegg_gene_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

ORA_KEGG_up_genes <- ORA_enrich(LANA_Up_genes, TERM2genes = m2tg)
dotplot(ORA_KEGG_up_genes)
ggsave(filename = "../BJAB_RNAseq/LANA_up_genes_BP.tiff")
ORA_KEGG_down_genes <- ORA_enrich(LANA_Down_genes, TERM2genes = m2tg)
dotplot(ORA_KEGG_down_genes)
ggsave(filename = "../BJAB_RNAseq/LANA_down_genes_BP.tiff")


ranked_gl <- step4_annotate_output_LANA_GFP$logFC
names(ranked_gl) <- step4_annotate_output_LANA_GFP$ensembl
c5BP_gene_sets = msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
msigdbr_list = split(x = c5BP_gene_sets$ensembl_gene, f = c5BP_gene_sets$gs_name)
c5BP_en <- GSEA_enrich(msigdbr_list, ranked_genes = ranked_gl)

topPathwaysUp <- c5BP_en[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- c5BP_en[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(msigdbr_list[topPathways], ranked_gl, c5BP_en,
              gseaParam=0.5)


# rank gene expression ---------------------------------------------------------------

gene_mat <- read.csv(file = "BJABs_LANA_OE_cpm.csv")
colnames(gene_mat)[1] <- "ensembl"

g_mat <- gene_mat[, 2:5]
rownames(g_mat) <- gene_mat$ensembl

boxplot(log2(g_mat))

gene_df <- g_mat %>% data_frame() %>% pivot_longer(names_to = "treatment", values_to = "count_per_m",
                                        cols = c("GFP1", "GFP2", "LANA1", "LANA2"))
ggplot(data = gene_df, aes(x=log2(count_per_m))) +
  geom_histogram(data=subset(gene_df, treatment == 'GFP1'), fill = "red", alpha = 0.2, bins = 2000) +
  geom_histogram(data=subset(gene_df, treatment == 'GFP2'), fill = "red", alpha = 0.2, bins = 2000) +
  geom_histogram(data=subset(gene_df, treatment == 'LANA1'), fill = "yellow", alpha = 0.2, bins = 2000) +
  geom_histogram(data=subset(gene_df, treatment == 'LANA2'), fill = "yellow", alpha = 0.2, bins = 2000)

g_mat %>% mutate(GFP = (GFP1 + GFP2)/2, LANA = (LANA1 + LANA2)/2) %>% ggplot(aes(x = log2(GFP), y = log2(LANA))) + geom_point()
hist(log2(g_mat[,1]), breaks = 1000)
hist(log2(g_mat[,2]), breaks = 1000)
hist(log2(g_mat[,3]), breaks = 1000)
hist(log2(g_mat[,4]), breaks = 1000)


