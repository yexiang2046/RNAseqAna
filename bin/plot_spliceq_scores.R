#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript plot_spliceq_scores.R <introns_with_peaks.bed> <introns_without_peaks.bed> [x_min] [x_max]")
}

# Load required libraries
if (!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
}
if (!require(ggridges)) {
    install.packages("ggridges")
    library(ggridges)
}

# Read input files
introns_with_peaks <- read.table(args[1], header=FALSE, 
                               col.names=c("chr", "start", "end", "intron_id", "score", 
                                         "strand", "spliceq_score", "peak_score", "peak_strand"))

introns_without_peaks <- read.table(args[2], header=FALSE,
                                  col.names=c("chr", "start", "end", "intron_id", "score", 
                                            "strand", "spliceq_score"))

# Add group labels
introns_with_peaks$group <- "With Peaks"
introns_without_peaks$group <- "Without Peaks"

# Combine data
combined_data <- rbind(
    data.frame(spliceq_score=introns_with_peaks$spliceq_score, group=introns_with_peaks$group),
    data.frame(spliceq_score=introns_without_peaks$spliceq_score, group=introns_without_peaks$group)
)

# Get counts for titles
with_peaks_count <- nrow(introns_with_peaks)
without_peaks_count <- nrow(introns_without_peaks)

# Get file names for titles
with_peaks_file <- basename(args[1])
without_peaks_file <- basename(args[2])

# Create violin plot
p1 <- ggplot(combined_data, aes(x=group, y=spliceq_score, fill=group)) +
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.2, alpha=0.7) +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    labs(
        title = paste("SPLICE-Q Score Distribution by Peak Overlap\n",
                     "With Peaks (", with_peaks_count, "): ", with_peaks_file, "\n",
                     "Without Peaks (", without_peaks_count, "): ", without_peaks_file),
        y = "SPLICE-Q Score"
    )

# Create ridge plot
p2 <- ggplot(combined_data, aes(x=spliceq_score, y=group, fill=group)) +
    geom_density_ridges(alpha=0.7) +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme_minimal() +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12, face="bold"),
        legend.position = "none"
    ) +
    labs(
        title = paste("SPLICE-Q Score Density Distribution\n",
                     "With Peaks (", with_peaks_count, "): ", with_peaks_file, "\n",
                     "Without Peaks (", without_peaks_count, "): ", without_peaks_file),
        x = "SPLICE-Q Score"
    )

# Set x-axis range if provided
if (length(args) >= 4) {
    x_min <- as.numeric(args[3])
    x_max <- as.numeric(args[4])
    if (!is.na(x_min) && !is.na(x_max)) {
        p2 <- p2 + scale_x_continuous(limits = c(x_min, x_max))
    }
}

# Save plots
ggsave("spliceq_scores_violin.pdf", p1, width=10, height=8)
ggsave("spliceq_scores_ridge.pdf", p2, width=10, height=6)

# Print summary statistics
cat("\nSummary Statistics:\n")
cat("\nIntrons with peaks:\n")
print(summary(introns_with_peaks$spliceq_score))
cat("\nIntrons without peaks:\n")
print(summary(introns_without_peaks$spliceq_score))

# Perform statistical test
cat("\nStatistical Test (Wilcoxon rank sum test):\n")
test_result <- wilcox.test(spliceq_score ~ group, data=combined_data)
print(test_result) 