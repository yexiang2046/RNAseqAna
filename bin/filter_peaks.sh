#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i INPUT_PEAKS -o OUTPUT_PEAKS [-p PVALUE] [-s SIGNAL] [-w WIDTH] [-n]"
    echo "Filter MACS2 narrowPeak output based on p-value, signal value, and peak width"
    echo ""
    echo "Arguments:"
    echo "  -i INPUT_PEAKS   Input narrowPeak file (required)"
    echo "  -o OUTPUT_PEAKS  Output filtered narrowPeak file (required)"
    echo "  -p PVALUE       P-value threshold (default: 1e-5)"
    echo "  -s SIGNAL       Signal value threshold (default: 10)"
    echo "  -w WIDTH        Peak width threshold in bp (default: 0, no filtering)"
    echo "  -n             Keep only peaks with negative p-values (default: false)"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hi:o:p:s:w:n" opt; do
    case $opt in
        h) usage ;;
        i) INPUT_PEAKS="$OPTARG" ;;
        o) OUTPUT_PEAKS="$OPTARG" ;;
        p) PVALUE="$OPTARG" ;;
        s) SIGNAL="$OPTARG" ;;
        w) WIDTH="$OPTARG" ;;
        n) NEGATIVE_ONLY=true ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$INPUT_PEAKS" ] || [ -z "$OUTPUT_PEAKS" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input file exists
if [ ! -f "$INPUT_PEAKS" ]; then
    echo "Error: Input peaks file does not exist: $INPUT_PEAKS"
    exit 1
fi

# Set default values
PVALUE=${PVALUE:-"1e-5"}
SIGNAL=${SIGNAL:-"10"}
WIDTH=${WIDTH:-"0"}

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_PEAKS")"

# Create temporary directory for intermediate files
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Function to check and fix tab formatting
check_and_fix_tabs() {
    local input_file="$1"
    local output_file="$2"
    local expected_fields="$3"
    
    # Check for tab issues
    if grep -q $'\t\t' "$input_file" || grep -q $'\t$' "$input_file"; then
        echo "Warning: Found extra tabs in $input_file. Fixing..."
        # Remove extra tabs and trailing tabs
        sed 's/\t\t/\t/g; s/\t$//' "$input_file" > "$output_file"
    else
        cp "$input_file" "$output_file"
    fi
    
    # Verify field count
    local field_count=$(awk -F'\t' 'NF != '"$expected_fields"' {print NR":"NF} END {print "Total lines:"NR}' "$output_file")
    if [[ $field_count != *"Total lines:"* ]]; then
        echo "Error: Found lines with incorrect number of fields:"
        echo "$field_count"
        exit 1
    fi
}

# Create R script for filtering
echo "Creating R script for peak filtering..."
cat > "$TEMP_DIR/filter_peaks.R" << 'EOF'
# Filter MACS2 peaks based on p-value, signal value, and peak width
library(tidyverse)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript filter_peaks.R input_file output_file pvalue_threshold signal_threshold width_threshold negative_only")
}

input_file <- args[1]
output_file <- args[2]
pvalue_threshold <- as.numeric(args[3])
signal_threshold <- as.numeric(args[4])
width_threshold <- as.numeric(args[5])
negative_only <- as.logical(args[6])

# Read the peaks file
peaks <- read.table(input_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Verify number of columns
if (ncol(peaks) != 10) {
  stop(sprintf("Error: Expected 10 columns, found %d columns", ncol(peaks)))
}

# Rename columns based on MACS2 output format
names(peaks) <- c("chr", "start", "end", "name", "score", "strand", 
                 "signalValue", "neglog10pvalue", "qValue", "peak")

# Calculate peak widths
peaks$width <- peaks$end - peaks$start

# Convert p-value threshold to -log10 scale
neglog10_threshold <- -log10(pvalue_threshold)

# Filter peaks based on criteria
filtered_peaks <- peaks %>%
  filter(signalValue >= signal_threshold)

if (width_threshold > 0) {
  filtered_peaks <- filtered_peaks %>%
    filter(width <= width_threshold)
}

if (negative_only) {
  filtered_peaks <- filtered_peaks %>%
    filter(neglog10pvalue >= neglog10_threshold)  # Higher -log10(p) means more significant
} else {
  filtered_peaks <- filtered_peaks %>%
    filter(abs(neglog10pvalue) >= neglog10_threshold)
}

# Write filtered peaks (full format)
write.table(filtered_peaks[, 1:10],  # Only write original columns
            file = output_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Write filtered peaks (BED format)
write.table(filtered_peaks[, 1:6],  # Write only BED columns
            file = paste0(output_file, ".bed"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Create output directory for plots
plot_dir <- dirname(output_file)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Create peak length distribution plot
pdf(file.path(plot_dir, "peak_length_distribution.pdf"), width = 10, height = 6)
ggplot(filtered_peaks, aes(x = width)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(x = "Peak Length (bp)",
       y = "Count",
       title = "Distribution of Peak Lengths",
       subtitle = sprintf("Mean: %.1f bp, Median: %.1f bp",
                        mean(filtered_peaks$width),
                        median(filtered_peaks$width)))
dev.off()

# Create signal value distribution plot
pdf(file.path(plot_dir, "signal_value_distribution.pdf"), width = 10, height = 6)
ggplot(filtered_peaks, aes(x = signalValue)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(x = "Signal Value",
       y = "Count",
       title = "Distribution of Peak Signal Values",
       subtitle = sprintf("Mean: %.2f, Median: %.2f",
                        mean(filtered_peaks$signalValue),
                        median(filtered_peaks$signalValue)))
dev.off()

# Calculate chromosome distribution
chr_stats <- filtered_peaks %>%
  group_by(chr) %>%
  summarise(
    peak_count = n(),
    percentage = n() / nrow(filtered_peaks) * 100
  ) %>%
  arrange(desc(peak_count))

# Print summary
cat(sprintf("\nPeak Filtering Summary:\n"))
cat(sprintf("Total peaks in input: %d\n", nrow(peaks)))
cat(sprintf("Peaks passing signal threshold (>= %.2f): %d\n", 
            signal_threshold, 
            sum(peaks$signalValue >= signal_threshold)))
if (width_threshold > 0) {
  cat(sprintf("Peaks passing width threshold (<= %d bp): %d\n", 
              width_threshold, 
              sum(peaks$width <= width_threshold)))
}
cat(sprintf("Peaks passing p-value threshold (<= %.2e): %d\n", 
            pvalue_threshold, 
            if(negative_only) {
              sum(peaks$neglog10pvalue >= neglog10_threshold)
            } else {
              sum(abs(peaks$neglog10pvalue) >= neglog10_threshold)
            }))
cat(sprintf("Final number of peaks: %d\n", nrow(filtered_peaks)))

# Print statistics
cat(sprintf("\nPeak Statistics:\n"))
cat(sprintf("Width Distribution:\n"))
cat(sprintf("  Min width: %d bp\n", min(peaks$width)))
cat(sprintf("  Max width: %d bp\n", max(peaks$width)))
cat(sprintf("  Mean width: %.1f bp\n", mean(peaks$width)))
cat(sprintf("  Median width: %.1f bp\n", median(peaks$width)))

cat(sprintf("\nP-value Distribution:\n"))
cat(sprintf("  Min -log10(p-value): %.2f\n", min(peaks$neglog10pvalue)))
cat(sprintf("  Max -log10(p-value): %.2f\n", max(peaks$neglog10pvalue)))
cat(sprintf("  Mean -log10(p-value): %.2f\n", mean(peaks$neglog10pvalue)))
cat(sprintf("  Median -log10(p-value): %.2f\n", median(peaks$neglog10pvalue)))
EOF

# Check and fix input files
echo "Checking input file formatting..."
check_and_fix_tabs "$INPUT_PEAKS" "$TEMP_DIR/input_peaks.bed" 10

# Run R script for filtering
echo "Running peak filtering..."
Rscript "$TEMP_DIR/filter_peaks.R" \
    "$TEMP_DIR/input_peaks.bed" \
    "$TEMP_DIR/filtered_peaks.bed" \
    "$PVALUE" \
    "$SIGNAL" \
    "$WIDTH" \
    "${NEGATIVE_ONLY:-false}"

# Check if filtering completed successfully
if [ $? -ne 0 ]; then
    echo "Error: Peak filtering failed"
    exit 1
fi

# Copy filtered peaks to output directory
cp "$TEMP_DIR/filtered_peaks.bed" "$OUTPUT_PEAKS"
cp "$TEMP_DIR/filtered_peaks.bed.bed" "$(dirname "$OUTPUT_PEAKS")/$(basename "$OUTPUT_PEAKS" .narrowPeak).bed"

# Print results
echo -e "\nPeak filtering complete!"
echo "Results written to:"
echo "- $OUTPUT_PEAKS (narrowPeak format)"
echo "- $(dirname "$OUTPUT_PEAKS")/$(basename "$OUTPUT_PEAKS" .narrowPeak).bed (BED format)"
echo "- $(dirname "$OUTPUT_PEAKS")/peak_length_distribution.pdf"
echo "- $(dirname "$OUTPUT_PEAKS")/signal_value_distribution.pdf"
echo "- $(dirname "$OUTPUT_PEAKS")/chromosome_distribution.pdf"

# Display summary if available
if [ -f "$OUTPUT_PEAKS" ]; then
    echo -e "\nFiltered peaks summary:"
    wc -l "$OUTPUT_PEAKS" | awk '{print "Number of peaks: " $1}'
fi 