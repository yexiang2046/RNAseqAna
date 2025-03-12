#!/bin/bash

# Help message
usage() {
    echo "Usage: $0 [-h] -p PEAKS_BED -r RMSK_BED -o OUTPUT_DIR"
    echo "Annotate peaks with RepeatMasker elements using bedtools annotate"
    echo ""
    echo "Arguments:"
    echo "  -p PEAKS_BED    Input peaks BED file (from Piranha)"
    echo "  -r RMSK_BED     RepeatMasker BED file (all_rmsk_hg38.bed)"
    echo "  -o OUTPUT_DIR   Output directory"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hp:r:o:" opt; do
    case $opt in
        h) usage ;;
        p) PEAKS="$OPTARG" ;;
        r) RMSK="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$PEAKS" ] || [ -z "$RMSK" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input files exist
if [ ! -f "$PEAKS" ]; then
    echo "Error: Peaks file does not exist: $PEAKS"
    exit 1
fi

if [ ! -f "$RMSK" ]; then
    echo "Error: RepeatMasker file does not exist: $RMSK"
    exit 1
fi

# Check if bedtools is installed
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed"
    exit 1
fi

# Extract sample ID from peaks filename
SAMPLE_ID=$(basename "$PEAKS" .bed | sed 's/_peaks$//')
echo "Processing sample: $SAMPLE_ID"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run annotateBed to get repeat element counts
echo "Running annotateBed for repeat element annotation..."
bedtools annotate -counts \
    -i "$PEAKS" \
    -files "$RMSK" \
    > "$OUTPUT_DIR/${SAMPLE_ID}_rmsk_counts.bed"

# Create R script for analysis
echo "Creating R script for analysis..."
cat > "$OUTPUT_DIR/analyze_rmsk.R" << 'EOF'
# Read and analyze RepeatMasker annotations
library(tidyverse)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript analyze_rmsk.R rmsk_counts.bed sample_id")
}

input_file <- args[1]
sample_id <- args[2]

# Read the counts file
data <- read.table(input_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Name columns appropriately
names(data) <- c("chr", "start", "end", "name", "score", "strand", "rmsk_count")

# Calculate summary statistics
summary_stats <- data %>%
  summarise(
    total_peaks = n(),
    peaks_with_repeats = sum(rmsk_count > 0),
    percent_with_repeats = mean(rmsk_count > 0) * 100,
    avg_repeats_per_peak = mean(rmsk_count),
    max_repeats = max(rmsk_count)
  )

# Create count distribution
count_dist <- data %>%
  count(rmsk_count) %>%
  mutate(
    percent = n / sum(n) * 100,
    rmsk_count = factor(rmsk_count, levels = as.character(sort(unique(rmsk_count))))
  )

# Write summary statistics
write.csv(summary_stats,
          file = file.path(dirname(input_file),
                          paste0(sample_id, "_rmsk_summary.csv")),
          row.names = FALSE)

# Create visualization
pdf(file.path(dirname(input_file),
              paste0(sample_id, "_rmsk_distribution.pdf")),
    width = 12, height = 8)

# Plot distribution of repeat counts
ggplot(count_dist, aes(x = rmsk_count, y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", percent)),
            vjust = -0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Number of Overlapping Repeat Elements",
       y = "Percentage of Peaks",
       title = paste("Distribution of Repeat Element Overlaps -", sample_id),
       subtitle = sprintf("Total peaks: %d, Peaks with repeats: %d (%.1f%%)",
                        summary_stats$total_peaks,
                        summary_stats$peaks_with_repeats,
                        summary_stats$percent_with_repeats))

dev.off()

# Print summary
cat(sprintf("\nRepeat Element Summary for %s:\n", sample_id))
print(summary_stats)
cat("\nDistribution of repeat counts:\n")
print(head(count_dist, 10))
EOF

# Run R script
echo "Running R analysis..."
Rscript "$OUTPUT_DIR/analyze_rmsk.R" \
    "$OUTPUT_DIR/${SAMPLE_ID}_rmsk_counts.bed" \
    "$SAMPLE_ID"

# Print results
echo -e "\nAnnotation complete!"
echo "Results written to:"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_rmsk_counts.bed"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_rmsk_summary.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_rmsk_distribution.pdf"

# Display summary if available
if [ -f "$OUTPUT_DIR/${SAMPLE_ID}_rmsk_summary.csv" ]; then
    echo -e "\nRepeat element summary:"
    cat "$OUTPUT_DIR/${SAMPLE_ID}_rmsk_summary.csv"
fi 