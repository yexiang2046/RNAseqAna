#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    if (!require(ChIPseeker)) BiocManager::install("ChIPseeker")
    if (!require(TxDb.Hsapiens.UCSC.hg38.knownGene)) BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
    if (!require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")
    if (!require(clusterProfiler)) BiocManager::install("clusterProfiler")
    if (!require(ggplot2)) install.packages("ggplot2")
    if (!require(tidyverse)) install.packages("tidyverse")
    
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(ggplot2)
    library(tidyverse)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript annotate_peaks.R <peaks_file> <output_dir> <sample_name>")
}

peaks_file <- args[1]
output_dir <- args[2]
sample_name <- args[3]

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read peaks file
peaks <- readPeakFile(peaks_file)

# Annotate peaks
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak_anno <- annotatePeak(peaks, 
                         TxDb = txdb,
                         tssRegion = c(-3000, 3000),
                         annoDb = "org.Hs.eg.db")

# Convert annotation to data frame
anno_df <- as.data.frame(peak_anno)

# Save annotated peaks
write.table(anno_df, 
            file = file.path(output_dir, paste0(sample_name, "_annotated_peaks.txt")),
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# Create distribution plots
pdf(file.path(output_dir, paste0(sample_name, "_peak_distribution.pdf")), 
    width = 10, height = 8)

# Distribution of peaks relative to TSS
plotDistToTSS(peak_anno, 
              title = paste("Distribution of Peaks Relative to TSS -", sample_name))

# Distribution of peaks by genomic features
plotAnnoPie(peak_anno, 
            title = paste("Distribution of Peaks by Genomic Features -", sample_name))

# Distribution of peaks by chromosome
plotAnnoBar(peak_anno, 
            title = paste("Distribution of Peaks by Chromosome -", sample_name))

dev.off()

# Create summary statistics
summary_stats <- list(
    total_peaks = length(peaks),
    peaks_in_promoter = sum(anno_df$annotation == "Promoter"),
    peaks_in_exon = sum(grepl("Exon", anno_df$annotation)),
    peaks_in_intron = sum(grepl("Intron", anno_df$annotation)),
    peaks_in_intergenic = sum(grepl("Intergenic", anno_df$annotation))
)

# Save summary statistics
write.table(as.data.frame(summary_stats),
            file = file.path(output_dir, paste0(sample_name, "_peak_stats.txt")),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Create gene annotation summary
gene_anno <- anno_df %>%
    filter(!is.na(SYMBOL)) %>%
    group_by(SYMBOL) %>%
    summarise(
        peak_count = n(),
        mean_score = mean(score),
        max_score = max(score),
        min_score = min(score)
    ) %>%
    arrange(desc(peak_count))

# Save gene annotation summary
write.table(gene_anno,
            file = file.path(output_dir, paste0(sample_name, "_gene_annotation.txt")),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Print summary
cat("\nPeak Annotation Summary:\n")
cat("------------------------\n")
cat(sprintf("Total peaks: %d\n", summary_stats$total_peaks))
cat(sprintf("Peaks in promoter regions: %d\n", summary_stats$peaks_in_promoter))
cat(sprintf("Peaks in exons: %d\n", summary_stats$peaks_in_exon))
cat(sprintf("Peaks in introns: %d\n", summary_stats$peaks_in_intron))
cat(sprintf("Peaks in intergenic regions: %d\n", summary_stats$peaks_in_intergenic))

# Create peak length distribution plot
pdf(file.path(output_dir, paste0(sample_name, "_peak_length_distribution.pdf")),
    width = 10, height = 6)

ggplot(anno_df, aes(x = width)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Peak Length (bp)",
         y = "Count",
         title = paste("Distribution of Peak Lengths -", sample_name),
         subtitle = sprintf("Mean: %.1f bp, Median: %.1f bp",
                          mean(anno_df$width),
                          median(anno_df$width)))

dev.off()

# Create peak score distribution plot
pdf(file.path(output_dir, paste0(sample_name, "_peak_score_distribution.pdf")),
    width = 10, height = 6)

ggplot(anno_df, aes(x = score)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(x = "Peak Score",
         y = "Count",
         title = paste("Distribution of Peak Scores -", sample_name),
         subtitle = sprintf("Mean: %.2f, Median: %.2f",
                          mean(anno_df$score),
                          median(anno_df$score)))

dev.off()

# Create chromosome distribution plot
chr_stats <- anno_df %>%
    group_by(seqnames) %>%
    summarise(
        peak_count = n(),
        percentage = n() / nrow(anno_df) * 100
    ) %>%
    arrange(desc(peak_count))

pdf(file.path(output_dir, paste0(sample_name, "_chromosome_distribution.pdf")),
    width = 12, height = 6)

ggplot(chr_stats, aes(x = reorder(seqnames, -peak_count), y = peak_count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%d (%.1f%%)", peak_count, percentage)),
              vjust = -0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Chromosome",
         y = "Number of Peaks",
         title = paste("Distribution of Peaks Across Chromosomes -", sample_name))

dev.off()

# Save chromosome statistics
write.table(chr_stats,
            file = file.path(output_dir, paste0(sample_name, "_chromosome_stats.txt")),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("\nResults written to:\n")
cat(sprintf("- %s\n", file.path(output_dir, paste0(sample_name, "_annotated_peaks.txt"))))
cat(sprintf("- %s\n", file.path(output_dir, paste0(sample_name, "_peak_stats.txt"))))
cat(sprintf("- %s\n", file.path(output_dir, paste0(sample_name, "_gene_annotation.txt"))))
cat(sprintf("- %s\n", file.path(output_dir, paste0(sample_name, "_peak_distribution.pdf"))))
cat(sprintf("- %s\n", file.path(output_dir, paste0(sample_name, "_peak_length_distribution.pdf"))))
cat(sprintf("- %s\n", file.path(output_dir, paste0(sample_name, "_peak_score_distribution.pdf"))))
cat(sprintf("- %s\n", file.path(output_dir, paste0(sample_name, "_chromosome_distribution.pdf"))))
cat(sprintf("- %s\n", file.path(output_dir, paste0(sample_name, "_chromosome_stats.txt")))) 