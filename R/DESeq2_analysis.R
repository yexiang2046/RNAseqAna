# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(RColorBrewer)

# Check for command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Insufficient arguments provided. Usage: Rscript this_script.R path/to/featureCounts_matrix.txt path/to/sample_info.txt path/to/output_folder/")
}

# Paths for data and output
countDataPath <- args[1]
sampleInfoPath <- args[2]
outputFolder <- args[3]

# Create the output folder if it does not exist
if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

# Load count data from featureCounts
countData <- read.table(countDataPath, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Load sample information
sampleInfo <- read.table(sampleInfoPath, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Check if column names of countData match row names of sampleInfo
if (!all(colnames(countData) == rownames(sampleInfo))) {
  stop("Mismatch between column names of count data and row names of sample info.")
}

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleInfo,
                              design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# PCA Plot
rld <- rlogTransformation(dds, blind = TRUE)  # Using rlog for normalization
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of rlog-transformed data") +
  theme_minimal()
ggsave(file.path(outputFolder, "PCA_plot.pdf"), plot = p)

# QC Plot: Dispersion Estimates
pdf(file.path(outputFolder, "Dispersion_plot.pdf"))
plotDispEsts(dds)
dev.off()

# Get results
res <- results(dds)

# Save results as an RDS file
saveRDS(res, file = file.path(outputFolder, "DESeq2_results.rds"))

# Print completion message
print("DESeq2 analysis completed. Results and plots saved in the specified output folder.")
