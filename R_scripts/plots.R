library(ggplot2)
library(pheatmap)




cpm_data <- read.csv("12630_RS_out/normalized_counts_CPM.csv")



files <- list.files("12630_RS_out/", 
                    pattern = "DEG_.*.csv",
                    include.dirs = TRUE,
                    full.names = TRUE)

th1_de_genes_file <- files[grep(".*Th1.*Th1.*", files)]
treg_de_genes_file <- files[grep(".*Treg.*Treg.*", files)]

th1_de_genes <- lapply(th1_de_genes_file, read.csv)
# combine the DEG genes pdfj < 0.01 and abs(logFC) > 1 into one dataframe
th1_de_genes <- do.call(rbind, th1_de_genes)
th1_de_genes <- th1_de_genes[th1_de_genes$FDR < 0.01 & abs(th1_de_genes$logFC) > 1, ]$gene_id
th1_de_genes <- unique(th1_de_genes)

treg_de_genes <- lapply(treg_de_genes_file, read.csv)
# combine the DEG genes pdfj < 0.01 into one dataframe
treg_de_genes <- do.call(rbind, treg_de_genes)
treg_de_genes <- treg_de_genes[treg_de_genes$FDR < 0.01 & abs(treg_de_genes$logFC) > 1, ]$gene_id
treg_de_genes <- unique(treg_de_genes)


cpm_data_th1_de <- cpm_data[cpm_data$gene_id %in% th1_de_genes, 3:22]
cpm_data_treg_de <- cpm_data[cpm_data$gene_id %in% treg_de_genes, 23:42]

metadata <- read.table("metaData_12630-RS.txt", sep = "\t", header = TRUE)
metadata$Temp <- as.factor(metadata$Temp)
colnames(cpm_data_th1_de) <- metadata$SampleId[1:20]
colnames(cpm_data_treg_de) <- metadata$SampleId[21:40]

annotation_col_th1 <- data.frame(
    cell_type = metadata$Celltype[1:20],
    temperature = metadata$Temp[1:20],
    genotype = metadata$genotype[1:20]
)
rownames(annotation_col_th1) <- metadata$SampleId[1:20]

annotation_col_treg <- data.frame(
    cell_type = metadata$Celltype[21:40],
    temperature = metadata$Temp[21:40],
    genotype = metadata$genotype[21:40]
)
rownames(annotation_col_treg) <- metadata$SampleId[21:40]

pdf("th1_de_heatmap.pdf")
    pheatmap(cpm_data_th1_de,
         scale = "row",
         cluster_rows = TRUE,
         annotation_col = annotation_col_th1,
         show_rownames = FALSE,
         cluster_cols = FALSE)
dev.off()

pdf("th1_heatmap_clustered.pdf")
    pheatmap(cpm_data[, 3:22],
         scale = "row",
         cluster_rows = FALSE,
         annotation_col = annotation_col_th1,
         show_rownames = FALSE,
         cluster_cols = FALSE)
dev.off()



pdf("treg_de_heatmap.pdf")
    pheatmap(cpm_data_treg_de,
         scale = "row",
         cluster_rows = FALSE,
         annotation_col = annotation_col_treg,
         show_rownames = FALSE,
         cluster_cols = FALSE)
dev.off()

pdf("treg_heatmap_clustered.pdf")
    pheatmap(cpm_data_treg_de,
         scale = "row",
         cluster_rows = TRUE,
         annotation_col = annotation_col_treg,
         show_rownames = FALSE,
         cluster_cols = FALSE)
dev.off()

# make volcano plots for all comparisons
dir.create("volcano_plots")
for (file in files) {
    de_genes <- read.csv(file)
    name <- gsub(".csv", "", basename(file))
    
    ggplot(de_genes, aes(x = logFC, y = -log10(FDR))) +
        geom_point() +
        geom_point(data = de_genes[de_genes$FDR < 0.01 & abs(de_genes$logFC) > 1, ], color = "red") +
        theme_minimal() +
        labs(title = name) +
        geom_vline(xintercept = 1, color = "red") +
        geom_vline(xintercept = -1, color = "red") +
        geom_hline(yintercept = 2, color = "blue") 
    ggsave(paste0("volcano_plots/", name, ".pdf"))
}
















