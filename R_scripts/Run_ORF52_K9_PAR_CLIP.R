
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
require(rtracklayer)


library(devtools)
load_all()
metaData <- data.frame(sampleID = c(paste0(7259, sep = '_', 1:4), paste0(7356, sep = '_', 1:4)),
                       Treatment = c("GFP_input", "GFP_input", "K9_input", "K9_input", 
                                     "K9_PARCLIP", "K9_PARCLIP", "ORF52_PARCLIP", "ORF52_PARCLIP"))






step1_get_counts_output <- step1_get_counts(featurecount_file = "PAR_CLIP_counts.txt",
                                            col_names = c("Geneid", "Chr", "Start", "End", "Strand", "Length", 
                                                          "GFP_input1", "GFP_input2", "K9_input1", "K9_input2", 
                                                          "K9_PARCLIP1", "K9_PARCLIP2", "ORF52_PARCLIP1", "ORF52_PARCLIP2"),
                                            counts_col = c("GFP_input1", "GFP_input2", "K9_input1", "K9_input2", 
                                                           "K9_PARCLIP1", "K9_PARCLIP2", "ORF52_PARCLIP1", "ORF52_PARCLIP2"))




group <- factor(metaData[, "Treatment"])
design <- model.matrix(~0+group)
con1 <- makeContrasts(groupK9_PARCLIP - groupK9_input, levels = design) # need to be generalized
con2 <- makeContrasts(groupORF52_PARCLIP - groupGFP_input, levels = design)

con <- list(con1, con2)
step2_edgeR_analysis_output <- step2_edgeR_analysis(counts = step1_get_counts_output, meta = metaData,
                                                    sample_column = "sampleID",
                                                    groups = "Treatment",
                                                    contrast_list = con,
                                                    logFC = 1, padj = 0.05)

step3_prep_annotation_output <- step3_prep_annotation(sp = "human",
                                                      IdList = rownames(step2_edgeR_analysis_output[[1]]$table),
                                                      info_types = c("entrez", "symbol", "ensembl"),
                                                      IdType = "ensembl")


step4_annotate_output_K9_PARCLIP <- step4_annotate(qlf = step2_edgeR_analysis_output[[1]],
                                               annotation = step3_prep_annotation_output,
                                               Id_col = "ensembl",
                                               file_name = "K9_CLIP_vs_GFP")
step4_annotate_output_52_PARCLIP <- step4_annotate(qlf = step2_edgeR_analysis_output[[2]],
                                                 annotation = step3_prep_annotation_output,
                                                 Id_col = "ensembl",
                                                 file_name = "ORF52_CLIP_vs_input")


hg38_annotation <- import.gff3(con = "GRCh38_HHV8.annotation.gff3")
hg38_tb <- data.frame(hg38_annotation)
hg38_tb$gene_id <- sapply(hg38_tb$gene_id, function(x){strsplit(x, "\\.")}[[1]][1])


K9_enriched_genes <- step4_annotate_output_K9_PARCLIP[step4_annotate_output_K9_PARCLIP$threshold == "Up",]
colnames(K9_enriched_genes)[1] <- "gene_id"
K9_enriched_genes <- left_join(K9_enriched_genes, hg38_tb, by = "gene_id")
K9_enriched_genes <- K9_enriched_genes[!duplicated(K9_enriched_genes$gene_id),]
K9_enriched_genes %>% write_csv(file = "K9_PARCLIP_enriched_genes.csv")

K9_type <- K9_enriched_genes$gene_type %>% table() %>% data.frame()
colnames(K9_type) <- c("gene_type", "counts")
K9_type %>% arrange(desc(counts)) %>% write_csv(file = "K9_PARCLIP_enriched_genes_type.csv")
K9_type %>% arrange(desc(counts)) %>%
  mutate(prop = counts / sum(K9_type$counts) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  ggplot(aes(x = "", y = counts, fill = gene_type)) +
  geom_bar(stat = "identity", width = 1) + coord_polar("y", start=0) +
  theme_void()
ggsave(filename = "K9_PARCLIP_enriched_genes_type.tiff")

ORF52_enriched_genes <- step4_annotate_output_52_PARCLIP[step4_annotate_output_52_PARCLIP$threshold == "Up",]
colnames(ORF52_enriched_genes)[1] <- "gene_id"
ORF52_enriched_genes <- left_join(ORF52_enriched_genes, hg38_tb, by = "gene_id")
ORF52_enriched_genes <- ORF52_enriched_genes[!duplicated(ORF52_enriched_genes$gene_id),]
ORF52_enriched_genes %>% write_csv(file = "ORF52_PARCLIP_enriched_genes.csv")
ORF52_type <- ORF52_enriched_genes$gene_type %>% table() %>% data.frame()
colnames(ORF52_type) <- c("gene_type", "counts")
ORF52_type %>% arrange(desc(counts)) %>% write_csv(file = "ORF52_PARCLIP_enriched_genes_type.csv")
ORF52_type %>% arrange(desc(counts)) %>%
  mutate(prop = counts / sum(ORF52_type$counts) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  ggplot(aes(x = "", y = counts, fill = gene_type)) +
  geom_bar(stat = "identity", width = 1) + coord_polar("y", start=0) +
  theme_void()
ggsave(filename = "ORF52_PARCLIP_enriched_genes_type.tiff")

kegg_gene_sets = msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
m2tg <- kegg_gene_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

ORA_KEGG_K9_genes <- ORA_enrich(K9_enriched_genes, TERM2genes = m2tg)
dotplot(ORA_KEGG_K9_genes)
ggsave(filename = "K9_PARCLIP_enriched_genes_BP.tiff")
ORA_KEGG_ORF52_genes <- ORA_enrich(ORF52_enriched_genes, TERM2genes = m2tg)
dotplot(ORA_KEGG_ORF52_genes)
ggsave(filename = "ORF52_PARCLIP_enriched_genes_BP.tiff")


