#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -t TREATMENT_BAM -c CONTROL_BAM -o OUTPUT_DIR -g GENOME_SIZE [-n NAME] [-q QVALUE] [-f FORMAT] [-B] [-C]"
    echo "Call peaks using MACS2"
    echo ""
    echo "Arguments:"
    echo "  -t TREATMENT_BAM  Treatment BAM file (required)"
    echo "  -c CONTROL_BAM    Control BAM file (required)"
    echo "  -o OUTPUT_DIR     Output directory"
    echo "  -g GENOME_SIZE    Genome size (e.g., hs for human, mm for mouse)"
    echo "  -n NAME          Name prefix for output files (default: macs2)"
    echo "  -q QVALUE        Q-value threshold (default: 0.05)"
    echo "  -f FORMAT        Output format (default: AUTO)"
    echo "  -B              Save fragment pileup signal"
    echo "  -C              Save control lambda signal"
    echo "  -h              Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "ht:c:o:g:n:q:f:BC" opt; do
    case $opt in
        h) usage ;;
        t) TREATMENT_BAM="$OPTARG" ;;
        c) CONTROL_BAM="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        g) GENOME_SIZE="$OPTARG" ;;
        n) NAME="$OPTARG" ;;
        q) QVALUE="$OPTARG" ;;
        f) FORMAT="$OPTARG" ;;
        B) SAVE_PILEUP=true ;;
        C) SAVE_LAMBDA=true ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$TREATMENT_BAM" ] || [ -z "$CONTROL_BAM" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$GENOME_SIZE" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input files exist
if [ ! -f "$TREATMENT_BAM" ]; then
    echo "Error: Treatment BAM file does not exist: $TREATMENT_BAM"
    exit 1
fi

if [ ! -f "$CONTROL_BAM" ]; then
    echo "Error: Control BAM file does not exist: $CONTROL_BAM"
    exit 1
fi

# Set default values
NAME=${NAME:-"macs2"}
QVALUE=${QVALUE:-"0.05"}
FORMAT=${FORMAT:-"BED"}

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Construct MACS2 command
MACS2_CMD="macs2 callpeak -t $TREATMENT_BAM -c $CONTROL_BAM -g $GENOME_SIZE -n $NAME --outdir $OUTPUT_DIR -q $QVALUE"

# Add optional parameters
if [ "$SAVE_PILEUP" = true ]; then
    MACS2_CMD="$MACS2_CMD -B"
fi

if [ "$SAVE_LAMBDA" = true ]; then
    MACS2_CMD="$MACS2_CMD -C"
fi

# Run MACS2
echo "Running MACS2 peak calling..."
echo "Command: $MACS2_CMD"
$MACS2_CMD

# Check if MACS2 completed successfully
if [ $? -ne 0 ]; then
    echo "Error: MACS2 peak calling failed"
    exit 1
fi

# Create R script for analysis
echo "Creating R script for peak analysis..."
cat > "$OUTPUT_DIR/analyze_peaks.R" << 'EOF'
# Read and analyze MACS2 peaks
library(tidyverse)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript analyze_peaks.R peaks_file sample_id")
}

peaks_file <- args[1]
sample_id <- args[2]

# Read the peaks file
peaks <- read.table(peaks_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Rename columns based on MACS2 output format
names(peaks) <- c("chr", "start", "end", "name", "score", "strand", 
                 "signalValue", "pValue", "qValue", "peak")

# Calculate peak lengths and statistics
peaks <- peaks %>%
  mutate(peak_length = end - start)

length_stats <- peaks %>%
  summarise(
    min_length = min(peak_length),
    max_length = max(peak_length),
    mean_length = mean(peak_length),
    median_length = median(peak_length),
    sd_length = sd(peak_length),
    total_peaks = n()
  )

# Write peak length statistics
write.csv(length_stats,
          file = file.path(dirname(peaks_file),
                          paste0(sample_id, "_peak_length_stats.csv")),
          row.names = FALSE)

# Create length distribution plot
pdf(file.path(dirname(peaks_file),
              paste0(sample_id, "_peak_length_distribution.pdf")),
    width = 10, height = 6)

ggplot(peaks, aes(x = peak_length)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(x = "Peak Length (bp)",
       y = "Count",
       title = paste("Distribution of Peak Lengths -", sample_id),
       subtitle = sprintf("Mean: %.1f bp, Median: %.1f bp",
                        length_stats$mean_length,
                        length_stats$median_length))

dev.off()

# Calculate signal value statistics
signal_stats <- peaks %>%
  summarise(
    min_signal = min(signalValue),
    max_signal = max(signalValue),
    mean_signal = mean(signalValue),
    median_signal = median(signalValue),
    sd_signal = sd(signalValue)
  )

# Write signal statistics
write.csv(signal_stats,
          file = file.path(dirname(peaks_file),
                          paste0(sample_id, "_signal_stats.csv")),
          row.names = FALSE)

# Create signal value distribution plot
pdf(file.path(dirname(peaks_file),
              paste0(sample_id, "_signal_distribution.pdf")),
    width = 10, height = 6)

ggplot(peaks, aes(x = signalValue)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(x = "Signal Value",
       y = "Count",
       title = paste("Distribution of Peak Signal Values -", sample_id),
       subtitle = sprintf("Mean: %.2f, Median: %.2f",
                        signal_stats$mean_signal,
                        signal_stats$median_signal))

dev.off()

# Calculate chromosome distribution
chr_stats <- peaks %>%
  group_by(chr) %>%
  summarise(
    peak_count = n(),
    percentage = n() / nrow(peaks) * 100
  ) %>%
  arrange(desc(peak_count))

# Write chromosome statistics
write.csv(chr_stats,
          file = file.path(dirname(peaks_file),
                          paste0(sample_id, "_chromosome_stats.csv")),
          row.names = FALSE)

# Create chromosome distribution plot
pdf(file.path(dirname(peaks_file),
              paste0(sample_id, "_chromosome_distribution.pdf")),
    width = 12, height = 6)

ggplot(chr_stats, aes(x = reorder(chr, -peak_count), y = peak_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%d (%.1f%%)", peak_count, percentage)),
            vjust = -0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Chromosome",
       y = "Number of Peaks",
       title = paste("Distribution of Peaks Across Chromosomes -", sample_id))

dev.off()

# Print summary
cat(sprintf("\nPeak Statistics for %s:\n", sample_id))
cat("\nLength Statistics:\n")
print(length_stats)

cat("\nSignal Value Statistics:\n")
print(signal_stats)

cat("\nChromosome Distribution:\n")
print(chr_stats)
EOF

# Run R analysis
echo "Running R analysis..."
Rscript "$OUTPUT_DIR/analyze_peaks.R" \
    "$OUTPUT_DIR/${NAME}_peaks.narrowPeak" \
    "$NAME"

# Print results
echo -e "\nPeak calling complete!"
echo "Results written to:"
echo "- $OUTPUT_DIR/${NAME}_peaks.narrowPeak"
echo "- $OUTPUT_DIR/${NAME}_peaks.xls"
echo "- $OUTPUT_DIR/${NAME}_summits.bed"
echo "- $OUTPUT_DIR/${NAME}_peak_length_stats.csv"
echo "- $OUTPUT_DIR/${NAME}_peak_length_distribution.pdf"
echo "- $OUTPUT_DIR/${NAME}_signal_stats.csv"
echo "- $OUTPUT_DIR/${NAME}_signal_distribution.pdf"
echo "- $OUTPUT_DIR/${NAME}_chromosome_stats.csv"
echo "- $OUTPUT_DIR/${NAME}_chromosome_distribution.pdf"

# Display summary if available
if [ -f "$OUTPUT_DIR/${NAME}_peak_length_stats.csv" ]; then
    echo -e "\nPeak length statistics:"
    cat "$OUTPUT_DIR/${NAME}_peak_length_stats.csv"
fi

if [ -f "$OUTPUT_DIR/${NAME}_signal_stats.csv" ]; then
    echo -e "\nSignal value statistics:"
    cat "$OUTPUT_DIR/${NAME}_signal_stats.csv"
fi 