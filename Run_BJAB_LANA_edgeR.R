
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

LANA_Up_genes <- step4_annotate_output_LANA_GFP[step4_annotate_output_LANA_GFP$threshold == "Up",]

LANA_Down_genes <- step4_annotate_output_LANA_GFP[step4_annotate_output_LANA_GFP$threshold == "Down",]

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
