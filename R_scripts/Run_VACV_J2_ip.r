

require(devtools)
require(readxl)
require(edgeR)
require(Organism.dplyr)
require(tidyverse)
load_all()
# library(Rsubread)


VACV_6h <- read_xlsx(path = "SAMPLE_SHEET/9814-YZ_NGS_UM_Library_Submission_Form_032723.xlsx", sheet = 2, skip = 3)
VACV_24h <- read_xlsx(path = "SAMPLE_SHEET/9221-YZ_NGS_UM_Library_Submission_Form__QC_form_Template.xlsx", sheet = 2, skip = 3)

colnames(VACV_6h)[1:5] <- c("SampleID", "Volume", "Conc", "Organism", "Description")
colnames(VACV_24h)[1:5] <- c("SampleID", "Volume", "Conc", "Organism", "Description")

metaData <- rbind(VACV_6h, VACV_24h)

metaData$Group <- c("VACV_input6h", "VACV_input6h", "E3Lmu_input6h", "E3Lmu_input6h", 
                    "VACV_IP6h", "VACV_IP6h", "E3Lmu_IP6h", "E3Lmu_IP6h",
                    "NonInf_input", "NonInf_input", "VACV_input24h", "VACV_input24h", "E3Lmu_input24h", "E3Lmu_input24h",
                    "NonInf_IP", "NonInf_IP", "VACV_IP24h", "VACV_IP24h", "E3Lmu_IP24h", "E3Lmu_IP24h")


VACV_human <- read.table(file = "TXT/count_genes_human.txt", header = TRUE)
VACV_virus <- read.table(file = "TXT/count_genes_VACV.txt", header = TRUE)

# used for viral gene analysis
cts <- rbind(VACV_human, VACV_virus)

# used for cellular gene analysis
# cts <- VACV_human


group <- factor(metaData$Group)

design <- model.matrix(~0 + group)
con1 <- makeContrasts(groupVACV_IP6h - groupVACV_input6h, levels = design)
con2 <- makeContrasts(groupE3Lmu_IP6h - groupE3Lmu_input6h, levels = design)
con3 <- makeContrasts(groupVACV_IP24h - groupVACV_input24h, levels = design)
con4 <- makeContrasts(groupE3Lmu_IP24h - groupE3Lmu_input24h, levels = design)
con5 <- makeContrasts(groupNonInf_IP - groupNonInf_input, levels = design)
#con <- list(con1, con2)

con6 <- makeContrasts(groupVACV_input6h - groupNonInf_input, levels = design)
con7 <- makeContrasts(groupE3Lmu_input6h - groupNonInf_input, levels = design)
con8 <- makeContrasts(groupE3Lmu_input6h - groupVACV_input6h, levels = design)

con9 <- makeContrasts(groupVACV_input24h - groupNonInf_input, levels = design)
con10 <- makeContrasts(groupE3Lmu_input24h - groupNonInf_input, levels = design)
con11 <- makeContrasts(groupE3Lmu_input24h - groupVACV_input24h, levels = design)

y<-DGEList(counts=cts, samples = metaData, group=group, annotation.columns = 1:6)
keep <- filterByExpr(y, group = group)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
y<-calcNormFactors(y)
y<-estimateGLMCommonDisp(y,design,verbose=TRUE)
y<-estimateGLMTrendedDisp(y,design)
y<-estimateGLMTagwiseDisp(y,design)
fit<-glmFit(y,design)


lrt1<-glmLRT(fit, contrast = con1)
lrt2<-glmLRT(fit, contrast = con2) 
lrt3<-glmLRT(fit, contrast = con3)
lrt4<-glmLRT(fit, contrast = con4)
lrt5<-glmLRT(fit, contrast = con5)

lrt6 <- glmLRT(fit, contrast = con6)
lrt7 <- glmLRT(fit, contrast = con7)
lrt8 <- glmLRT(fit, contrast = con8)

lrt9 <- glmLRT(fit, contrast = con9)
lrt10 <- glmLRT(fit, contrast = con10)
lrt11 <- glmLRT(fit, contrast = con11)


V6h_tb <- lrt1@.Data[[15]]
colnames(lrt1[[12]])
src <- src_ucsc(organism = "human")
IdList <- lrt1[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(V6h_tb) <- IdList
V6h_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
V6h_tb <- left_join(V6h_tb, sym, by = "ensembl", multiple = "first")
V6h_tb$padj <- p.adjust(V6h_tb$PValue, method = "BH")
V6h_tb %>% filter(logCPM > 0) %>% write.csv(file = "VACV_6h_J2IP_enriched_transcripts.csv")



library(ggplot2)
library(ggrepel)
V6h_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("WT6h-----------logFC-----------WT6h_IP") +
  ggtitle("VACV_WT6h_J2IP")

V6h_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_hline(yintercept = 2) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("WT6h-----------logFC-----------WT6h_IP") +
  ggtitle("VACV_WT6h_J2IP")

V6h_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% write.csv(file = "VACV_WT_viral_6h_J2_IP_viral_CDS.csv")

E3L6h_tb <- lrt2@.Data[[15]]
colnames(lrt2[[12]])
src <- src_ucsc(organism = "human")
IdList <- lrt2[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(E3L6h_tb) <- IdList
E3L6h_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
E3L6h_tb <- left_join(E3L6h_tb, sym, by = "ensembl", multiple = "first")
E3L6h_tb$padj <- p.adjust(E3L6h_tb$PValue, method = "BH")
E3L6h_tb %>% filter(logCPM > 0) %>% write.csv(file = "VACVdE3L_6h_J2IP_enriched_transcripts.csv")

E3L6h_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(J2 IP/dE3L 6h)") +
  ggtitle("VACV_dE3L_6h_J2IP")

E3L6h_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_hline(yintercept = 2) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("dE3L6h-----------logFC-----------dE3L6h_IP") +
  ggtitle("VACV_dE3L6h_J2IP")

E3L6h_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% write.csv(file = "VACV_dE3L_6h_J2_IP_viral_CDS.csv")

Uninf_tb <- lrt5@.Data[[15]]
colnames(lrt5[[12]])
src <- src_ucsc(organism = "human")
IdList <- lrt5[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(Uninf_tb) <- IdList
Uninf_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
Uninf_tb <- left_join(Uninf_tb, sym, by = "ensembl", multiple = "first")
Uninf_tb$padj <- p.adjust(Uninf_tb$PValue, method = "BH")
Uninf_tb %>% filter(logCPM > 0) %>% write.csv(file = "Uninfected_J2IP_enriched_transcripts.csv")

Uninf_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_hline(yintercept = 2) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("Mock-----------logFC-----------Mock_IP") +
  ggtitle("VACV_Mock_J2IP")

Uninf_tb_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% write.csv(file = "Mock_J2_IP_viral_CDS.csv")

Uninf_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(J2 IP/Uninf input)") +
  ggtitle("uninf_J2IP")

Uninf_tb_sig<- Uninf_tb %>% filter(logCPM > 0) %>% filter(logFC > 1 & padj < 0.01) %>% filter(str_detect(ensembl, "^cds", negate = TRUE))
Uninf_tb_sig %>% write.csv(file = "Uninfected_J2_IP_host_transcripts.csv")






VACV_24h_tb <- lrt3@.Data[[15]]
colnames(lrt3[[12]])
src <- src_ucsc(organism = "human")
IdList <- lrt3[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(VACV_24h_tb) <- IdList
VACV_24h_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
VACV_24h_tb <- left_join(VACV_24h_tb, sym, by = "ensembl", multiple = "first")
VACV_24h_tb$padj <- p.adjust(VACV_24h_tb$PValue, method = "BH")
VACV_24h_tb %>% filter(logCPM > 0) %>% write.csv(file = "VACV_24h_J2IP_enriched_transcripts.csv")



VACV_24h_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_hline(yintercept = 2) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("WT24h-----------logFC-----------WT24_IP") +
  ggtitle("VACV_WT_24h_J2IP")

VACV_24h_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% write.csv(file = "VACV_WT_24h_J2_IP_viral_CDS.csv")

VACV_24h_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(J2 IP/WT 24h)") +
  ggtitle("VACV_WT_24h_J2IP")


VACVdE3L_24h_tb <- lrt4@.Data[[15]]
src <- src_ucsc(organism = "human")
IdList <- lrt4[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(VACVdE3L_24h_tb) <- IdList
VACVdE3L_24h_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
VACVdE3L_24h_tb<- left_join(VACVdE3L_24h_tb, sym, by = "ensembl", multiple = "first")
VACVdE3L_24h_tb$padj <- p.adjust(VACVdE3L_24h_tb$PValue, method = "BH")
VACVdE3L_24h_tb%>% filter(logCPM > 0) %>% write.csv(file = "VACVdE3L_24h_J2IP_enriched_transcripts.csv")


VACVdE3L_24h_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_hline(yintercept = 2) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("dE3L24h-----------logFC-----------dE3L24_IP") +
  ggtitle("VACV_dE3L_24h_J2IP")

VACVdE3L_24h_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% 
  filter(logCPM > 0 ) %>% 
  write.csv(file = "VACV_dE3L_24h_J2_IP_viral_CDS.csv")

VACVdE3L_24h_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(logFC > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(J2 IP/dE3L 24h)") +
  ggtitle("VACV_dE3L_24h_J2IP")

VACV6h_uninf_tb <- lrt6@.Data[[15]]
src <- src_ucsc(organism = "human")
IdList <- lrt6[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(VACV6h_uninf_tb) <- IdList
VACV6h_uninf_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
VACV6h_uninf_tb <- left_join(VACV6h_uninf_tb, sym, by = "ensembl", multiple = "first")
VACV6h_uninf_tb$padj <- p.adjust(VACV6h_uninf_tb$PValue, method = "BH")
# VACV6h_uninf_tb%>% filter(logCPM > 0) %>% write.csv(file = "VACV6h_uninf_enriched_cellular_transcripts.csv")

VACV6h_uninf_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(abs(logFC) > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) + 
  geom_hline(yintercept = 2) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(WT6h/uninf)") +
  ggtitle("WT6h_DEGs")


VACV6h_uninf_tb_sig_up <- VACV6h_uninf_tb[VACV6h_uninf_tb$padj < 0.01 & VACV6h_uninf_tb$logFC > 1, ]
VACV6h_uninf_tb_sig_down <- VACV6h_uninf_tb[VACV6h_uninf_tb$padj < 0.01 & VACV6h_uninf_tb$logFC < -1, ]

ekg <- enrichGO(gene = VACV6h_uninf_tb[VACV6h_uninf_tb$padj < 0.01,]$entrez, OrgDb = org.Hs.eg.db, ont = "BP")
ekg <- simplify(ekg, cutoff = 0.7)
dotplot(ekg)



VACV6hdE3L_uninf_tb <- lrt7@.Data[[15]]
src <- src_ucsc(organism = "human")
IdList <- lrt7[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(VACV6hdE3L_uninf_tb) <- IdList
VACV6hdE3L_uninf_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
VACV6hdE3L_uninf_tb <- left_join(VACV6hdE3L_uninf_tb, sym, by = "ensembl", multiple = "first")
VACV6hdE3L_uninf_tb$padj <- p.adjust(VACV6hdE3L_uninf_tb$PValue, method = "BH")
# VACV6hdE3L_uninf_tb%>% filter(logCPM > 0) %>% write.csv(file = "VACV6hdE3L_uninf_enriched_cellular_transcripts.csv")

VACV6hdE3L_uninf_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(abs(logFC) > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) + 
  geom_hline(yintercept = 2) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(dE3L6h/uninf)") +
  ggtitle("dE3L6h_DEGs")


VACV6hdE3L_uninf_tb_sig_up <- VACV6hdE3L_uninf_tb[VACV6hdE3L_uninf_tb$padj < 0.01 & VACV6hdE3L_uninf_tb$logFC > 1,]
VACV6hdE3L_uninf_tb_sig_down <- VACV6hdE3L_uninf_tb[VACV6hdE3L_uninf_tb$padj < 0.01 & VACV6hdE3L_uninf_tb$logFC < -1,]

ekg <- enrichGO(gene = VACV6hdE3L_uninf_tb[VACV6hdE3L_uninf_tb$padj < 0.01,]$entrez, OrgDb = org.Hs.eg.db, ont = "BP")
ekg <- simplify(ekg, cutoff = 0.7)
dotplot(ekg)

# 6h scatter plot
VACV6_tb <- left_join(VACV6h_uninf_tb, VACV6hdE3L_uninf_tb, by = "ensembl", multiple = "first", suffix = c(".WT", ".dE3L"))
label_data <- VACV6_tb %>% subset(abs(logFC.WT) > 5) %>%
  subset(abs(logFC.dE3L) > 5) %>% subset(-log10(padj.WT) > 5) %>% subset(-log10(padj.dE3L) > 5)
  
VACV6_tb %>% filter(str_detect(ensembl, "^cds", negate = TRUE)) %>% 
  ggplot(aes(x = logFC.WT, y = logFC.dE3L, label = symbol.WT)) + geom_point() + theme_classic() + xlim(c(-15, 15)) +
  ylim(c(-15, 15)) + geom_text_repel(data = label_data, aes(x = logFC.WT, y = logFC.dE3L, label = symbol.WT),force = 5) + ggtitle("VACV_6h_induced_genes") + ylab("dE3L/Mock") + xlab("Parental/Mock") + geom_abline(intercept =0 , slope = 1, color = "red")

# table generation ----------

VACV6hdE3L_WT_tb <- lrt8@.Data[[15]]
src <- src_ucsc(organism = "human")
IdList <- lrt8[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(VACV6hdE3L_WT_tb) <- IdList
VACV6hdE3L_WT_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
VACV6hdE3L_WT_tb <- left_join(VACV6hdE3L_WT_tb, sym, by = "ensembl", multiple = "first")
VACV6hdE3L_WT_tb$padj <- p.adjust(VACV6hdE3L_WT_tb$PValue, method = "BH")
# VACV6hdE3L_WT_tb%>% filter(logCPM > 0) %>% write.csv(file = "cellular_VACV6hdE3L_WT_enriched_transcripts.csv")
VACV6hdE3L_WT_tb%>% filter(logCPM > 0) %>% write.csv(file = "VACV6hdE3L_WT_enriched_viral_transcripts.csv")



label_data <- VACV6hdE3L_WT_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE)) 

label_data %>% write.csv(file = "VACV6hh_dE3L_WT_viral_CDS.csv")


VACV6hdE3L_WT_tb %>% filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% ggplot(aes(x = logFC, y = -log10(padj), label = ensembl, colour = ifelse(abs(logFC) > 1 & padj < 0.01, "DEGS", "NS"))) + 
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) + 
  geom_hline(yintercept = 2) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(dE3L6h/WT)") +
  ggtitle("dE3L6h_WT_DEGs")  
  # + geom_text_repel(data = label_data, aes(x = logFC, y = -log10(padj), label = ensembl)) + ggtitle("VACV_6h_induced_genes") + ylab("dE3L/Mock") + xlab("Parental/Mock")

ekg <- enrichGO(gene = VACV6hdE3L_WT_tb[VACV6hdE3L_WT_tb$padj < 0.01,]$entrez, OrgDb = org.Hs.eg.db, ont = "BP")
ekg <- simplify(ekg, cutoff = 0.7)
dotplot(ekg)


# 24h scatter plot -----------------
VACV24h_uninf_tb <- lrt9@.Data[[15]]
src <- src_ucsc(organism = "human")
IdList <- lrt9[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(VACV24h_uninf_tb) <- IdList
VACV24h_uninf_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
VACV24h_uninf_tb <- left_join(VACV24h_uninf_tb, sym, by = "ensembl", multiple = "first")
VACV24h_uninf_tb$padj <- p.adjust(VACV24h_uninf_tb$PValue, method = "BH")
# VACV24h_uninf_tb%>% filter(logCPM > 0) %>% write.csv(file = "VACV24h_uninf_enriched_cellular_transcripts.csv")

VACV24h_uninf_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(abs(logFC) > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) + 
  geom_hline(yintercept = 2) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(WT24h/uninf)") +
  ggtitle("WT24h_DEGs")

VACV24h_uninf_tb_sig_up <- VACV24h_uninf_tb[VACV24h_uninf_tb$padj < 0.01 & VACV24h_uninf_tb$logFC > 1,]
VACV24h_uninf_tb_sig_down <- VACV24h_uninf_tb[VACV24h_uninf_tb$padj < 0.01 & VACV24h_uninf_tb$logFC < -1,]


ekg <- enrichGO(gene = VACV24h_uninf_tb[VACV24h_uninf_tb$padj < 0.01,]$entrez, OrgDb = org.Hs.eg.db, ont = "BP")
ekg <- simplify(ekg, cutoff = 0.7)
dotplot(ekg)



VACV24hdE3L_uninf_tb <- lrt10@.Data[[15]]
src <- src_ucsc(organism = "human")
IdList <- lrt10[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(VACV24hdE3L_uninf_tb) <- IdList
VACV24hdE3L_uninf_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
VACV24hdE3L_uninf_tb <- left_join(VACV24hdE3L_uninf_tb, sym, by = "ensembl", multiple = "first")
VACV24hdE3L_uninf_tb$padj <- p.adjust(VACV24hdE3L_uninf_tb$PValue, method = "BH")
VACV24hdE3L_uninf_tb%>% filter(logCPM > 0) %>% write.csv(file = "VACV24hdE3L_uninf_enriched_cellular_transcripts.csv")

VACV24hdE3L_uninf_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(abs(logFC) > 1 & PValue < 0.01, "DEGS", "NS"))) + geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) + 
  geom_hline(yintercept = 2) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(dE3L24h/uninf)") +
  ggtitle("dE3L24h_DEGs")

VACV24hdE3L_uninf_tb_sig_up <- VACV24hdE3L_uninf_tb[VACV24hdE3L_uninf_tb$padj < 0.01 & VACV24hdE3L_uninf_tb$logFC > 1,]
VACV24hdE3L_uninf_tb_sig_down <- VACV24hdE3L_uninf_tb[VACV24hdE3L_uninf_tb$padj < 0.01 & VACV24hdE3L_uninf_tb$logFC < -1,]


ekg <- enrichGO(gene = VACV24hdE3L_uninf_tb[VACV24hdE3L_uninf_tb$padj < 0.01,]$entrez, OrgDb = org.Hs.eg.db, ont = "BP")
ekg <- simplify(ekg, cutoff = 0.7)
dotplot(ekg)


VACV24hdE3L_WT_tb <- lrt11@.Data[[15]]
src <- src_ucsc(organism = "human")
IdList <- lrt10[[12]]$Geneid
IdList <- sapply(IdList, function(x){strsplit(x, "\\.")}[[1]][[1]])
rownames(VACV24hdE3L_WT_tb) <- IdList
VACV24hdE3L_WT_tb$ensembl <- IdList
sym <- select(src, keys = IdList, 
              columns = c("entrez", "symbol", "ensembl"),
              keytype = "ensembl")
colnames(sym) <- c("ensembl", "entrez", "symbol")
VACV24hdE3L_WT_tb <- left_join(VACV24hdE3L_WT_tb, sym, by = "ensembl", multiple = "first")
VACV24hdE3L_WT_tb$padj <- p.adjust(VACV24hdE3L_WT_tb$PValue, method = "BH")
# VACV24hdE3L_WT_tb%>% filter(logCPM > 0) %>% write.csv(file = "VACV24hdE3L_WT_enriched_cellular_transcripts.csv")
VACV24hdE3L_WT_tb_sig <- VACV24hdE3L_WT_tb[VACV24hdE3L_WT_tb$padj < 0.01,]

label_data <- VACV24hdE3L_WT_tb %>% 
  filter(str_detect(ensembl, "^cds", negate = FALSE))
label_data %>% write.csv(file = "VACV24h_dE3L_WT_viral_CDS.csv")

VACV24hdE3L_WT_tb %>% filter(str_detect(ensembl, "^cds", negate = FALSE)) %>% ggplot(aes(x = logFC, y = -log10(padj), label = ensembl, colour = ifelse(abs(logFC) > 1 & padj < 0.01, "DEGS", "NS"))) + 
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) + 
  geom_hline(yintercept = 2) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + 
  xlab("log(dE3L24h/WT)") +
  ggtitle("dE3L24h_WT_DEGs")  


# scatter plot
VACV24h_tb <- left_join(VACV24h_uninf_tb, VACV24hdE3L_uninf_tb, by = "ensembl", multiple = "first", suffix = c(".WT", ".dE3L"))

label_data <- VACV24h_tb %>% subset(abs(logFC.WT) > 5) %>%
  subset(abs(logFC.dE3L) > 5) %>% subset(-log10(padj.WT) > 5) %>% subset(-log10(padj.dE3L) > 5)
VACV24h_tb %>% filter(str_detect(ensembl, "^cds", negate = TRUE)) %>% 
  ggplot(aes(x = logFC.WT, y = logFC.dE3L, label = symbol.WT)) + geom_point() + theme_classic() + xlim(c(-15, 15)) +
  ylim(c(-15, 15)) + geom_text_repel(data = label_data, aes(x = logFC.WT, y = logFC.dE3L, label = symbol.WT),force = 5) + ggtitle("VACV_24h_induced_genes") + ylab("dE3L/Mock") + xlab("Parental/Mock") + geom_abline(intercept =0 , slope = 1, color = "red")



library(ggplot2)
library(ggrepel)

VACV6hdE3L_WT_tb_viral <- VACV6hdE3L_WT_tb %>% filter(grepl("^cds", ensembl))

VACV6hdE3L_WT_tb_viral %>% write.csv(file = "Differential_viral_transcripts_6h.csv")

VACV6hdE3L_WT_tb_viral %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = ensembl, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = VACV6hdE3L_WT_tb_viral, aes(x = logFC, y = -log10(padj), label = ensembl, max.overlaps = 200)) + xlim(c(-20, 20)) +
  xlab("VACV_WT-----------logFC6h-----------VACV_dE3L") +
  ggtitle("VACVdE3L_WT_6h_viral")

VACV24hdE3L_WT_tb_viral <- VACV24hdE3L_WT_tb %>% filter(grepl("^cds", ensembl))
VACV24hdE3L_WT_tb_viral %>% write.csv(file = "Differential_viral_transcripts_24h.csv")

VACV24hdE3L_WT_tb_viral %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = ensembl, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = VACV6hdE3L_WT_tb_viral, aes(x = logFC, y = -log10(padj), label = ensembl, max.overlaps = 200)) + xlim(c(-20, 20)) +
  xlab("VACV_WT-----------logFC24h-----------VACV_dE3L") +
  ggtitle("VACVdE3L_WT_24h_viral")

subdat <- VACV6h_uninf_tb %>% filter(logCPM > 0) %>% subset(-log10(padj) > 30)
VACV6h_uninf_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = subdat, aes(x = logFC, y = -log10(padj), label = symbol, max.overlaps = 200)) + xlim(c(-20, 20)) +
  xlab("Uninf-----------logFC-----------VACV_WT_6h") +
  ggtitle("VACV_WT6h_Uninf")

subdat <- VACV6hdE3L_uninf_tb %>% filter(logCPM > 0) %>% subset(-log10(padj) > 30)
VACV6hdE3L_uninf_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = subdat, aes(x = logFC, y = -log10(padj), label = symbol, max.overlaps = 200)) + xlim(c(-20, 20)) +
  xlab("Uninf-----------logFC-----------VACV_dE3L_6h") +
  ggtitle("VACV_dE3L_6h_Uninf")

subdat <- VACV6hdE3L_WT_tb %>% subset(-log10(padj) > 10)
VACV6hdE3L_WT_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + 
  geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = subdat, aes(x = logFC, y = -log10(padj), label = ensembl, max.overlaps = 200)) + xlim(c(-20, 20)) +
  xlab("WT-----------logFC-----------VACV_dE3L_6h") +
  ggtitle("VACV_dE3L_6h_WT")


subdat <- V6h_tb %>% filter(logCPM > 0) %>% subset(-log10(PValue) > 80)
  V6h_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(PValue), label = symbol, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = subdat, aes(x = logFC, y = -log10(PValue), label = ensembl, max.overlaps = 200)) + 
  xlab("-----------logFC-----------IPed") +
  ggtitle("VACV_6h_J2IP")


  subdat <- E3L6h_tb %>% filter(logCPM > 0) %>% subset(-log10(PValue) > 80)
  E3L6h_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(PValue), label = symbol, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = subdat, aes(x = logFC, y = -log10(PValue), label = symbol, max.overlaps = 200)) + 
  xlab("-----------logFC-----------IPed") +
  ggtitle("VACVdE3L_6h_J2IP")


  subdat <- Uninf_tb %>% filter(logCPM > 0) %>% subset(-log10(PValue) > 80)
  Uninf_tb %>% filter(logCPM >0) %>% ggplot(aes(x = logFC, y = -log10(PValue), label = symbol, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = subdat, aes(x = logFC, y = -log10(PValue), label = symbol, max.overlaps = 200)) + 
  xlab("-----------logFC-----------IPed") +
  ggtitle("Unifected_J2IP")


  subdat <- VACV_24h_tb %>% filter(logCPM > 0) %>% subset(VACV_24h_tb, -log10(PValue) > 80)
  VACV_24h_tb %>% filter(logCPM > 0) %>% ggplot(aes(x = logFC, y = -log10(PValue), label = symbol, colour = ifelse(abs(logFC) > 0.5 & PValue < 0.00001, "DEGS", "NS"))) + geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = subdat, aes(x = logFC, y = -log10(PValue), label = ensembl, max.overlaps = 200)) + 
  xlab("-----------logFC-----------IPed") +
  ggtitle("VACV_24h_J2IP")


  subdat <- VACVdE3L_24h_tb %>% filter(logCPM > 0) %>%
   subset(-log10(padj) > 50 & logFC > 0)
  VACVdE3L_24h_tb %>% filter(logCPM > 0)%>% ggplot(aes(x = logFC, y = -log10(padj), label = symbol, colour = ifelse(abs(logFC) > 0.5 & padj < 1e-30, "DEGS", "NS"))) + geom_point() + theme_classic() +
  scale_color_manual(values = c("darkred", "gray80")) +
  labs(fill = "Genes") + geom_text_repel(data = subdat, aes(x = logFC, y = -log10(padj), label = symbol, max.overlaps = 200)) + 
  xlab("-----------logFC-----------IPed") +
  ggtitle("VACVdE3L_24h_J2IP")



