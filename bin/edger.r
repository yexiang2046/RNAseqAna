#!/usr/bin/env Rscript

# Load necessary libraries
library(edgeR)
library(optparse)
library(ggplot2)
library(factoextra)

# Define command-line options
option_list <- list(
  make_option(c("-c", "--counts"), type = "character", help = "Path to counts file"),
  make_option(c("-m", "--metadata"), type = "character", help = "Path to metadata file"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory for results")
)

# Parse command-line options
opt <- parse_args(OptionParser(option_list = option_list))

# Load count data and metadata
counts <- read.delim(opt$counts, header = TRUE, skip = 1)
colnames(counts) <- sapply(colnames(counts), function(x){strsplit(x, "_")}[[1]][1])
colnames(counts) <- gsub("\\.", "-", colnames(counts))
colnames(counts) <- gsub("^X", "", colnames(counts))
# Extract gene expression counts (columns 7 onwards) 
counts_matrix <- counts[,7:ncol(counts)]
# Keep annotation columns separately if needed
gene_annotations <- counts[,1:6]
metaData <- read.table(opt$metadata, header = TRUE)

# order counts_matrix to metaData sample order with SampleId
counts_matrix <- counts_matrix[, metaData$SampleId]



# Create DGEList object
group <- factor(metaData$group)
y <- DGEList(counts = counts_matrix, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)

# Create design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit the model
fit <- glmQLFit(y, design)

# Perform pairwise comparisons
group_levels <- levels(group)
comparison_results <- list()

for (i in 1:(length(group_levels) - 1)) {
  for (j in (i + 1):length(group_levels)) {
    contrast_name <- paste(group_levels[i], "vs", group_levels[j], sep = "_")
    contrast <- makeContrasts(contrasts = paste(group_levels[i], "-", group_levels[j]), levels = design)
    qlf <- glmQLFTest(fit, contrast = contrast)
    results <- topTags(qlf, n = Inf)
    comparison_results[[contrast_name]] <- results
    write.csv(results, file = file.path(opt$output, paste0("DEG_", contrast_name, ".csv")))
  }
}

# Plot PCA
pca_data <- prcomp(t(cpm(y, log = TRUE)), scale. = TRUE)
pca_plot <- fviz_pca_ind(pca_data, 
                         label = "none", 
                         habillage = group, 
                         addEllipses = FALSE) +
  ggtitle("PCA of Samples") +
  theme_classic()
ggsave(filename = file.path(opt$output, "PCA_plot.png"), plot = pca_plot)

# Bar graph of differentially expressed genes for each comparison
for (contrast_name in names(comparison_results)) {
  de_genes <- comparison_results[[contrast_name]]$table
  de_genes$threshold <- as.factor(
    ifelse(de_genes$FDR < 0.05 & abs(de_genes$logFC) > 1,
           "Significant", 
           "Not Significant")
  )
  de_counts <- table(de_genes$threshold)
  bar_plot <- ggplot(as.data.frame(de_counts), 
                     aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = paste("Number of Differentially Expressed Genes:", 
                      contrast_name),
         x = "Group",
         y = "Count")
  ggsave(filename = file.path(opt$output, 
                             paste0("DEG_barplot_", contrast_name, ".png")),
         plot = bar_plot)
}