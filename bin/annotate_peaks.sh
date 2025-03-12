#!/bin/bash

# Help message
usage() {
    echo "Usage: $0 [-h] -p PEAKS_BED -f FEATURES_DIR -o OUTPUT_DIR"
    echo "Annotate peaks with genomic features using bedtools annotate"
    echo ""
    echo "Arguments:"
    echo "  -p PEAKS_BED    Input peaks BED file (from Piranha)"
    echo "  -f FEATURES_DIR Directory containing feature BED files"
    echo "  -o OUTPUT_DIR   Output directory"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hp:f:o:" opt; do
    case $opt in
        h) usage ;;
        p) PEAKS="$OPTARG" ;;
        f) FEATURES_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$PEAKS" ] || [ -z "$FEATURES_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input files exist
if [ ! -f "$PEAKS" ]; then
    echo "Error: Peaks file does not exist: $PEAKS"
    exit 1
fi

if [ ! -d "$FEATURES_DIR" ]; then
    echo "Error: Features directory does not exist: $FEATURES_DIR"
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

# List of feature types to process
FEATURES=("cds" "exons" "five_prime_utr" "genes" "introns" "three_prime_utr" "transcripts")

# Check if all feature files exist
for feature in "${FEATURES[@]}"; do
    if [ ! -f "$FEATURES_DIR/${feature}.bed" ]; then
        echo "Error: Feature file not found: $FEATURES_DIR/${feature}.bed"
        exit 1
    fi
done

# Create feature files list for annotateBed
FEATURE_FILES=""
for feature in "${FEATURES[@]}"; do
    FEATURE_FILES="$FEATURE_FILES $FEATURES_DIR/${feature}.bed"
done

# First, get gene annotations using bedtools intersect
echo "Getting gene annotations..."
bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/genes.bed" -wao | \
    awk 'BEGIN{OFS="\t"} {
        # Store gene info for each peak
        peak=$1"_"$2"_"$3;
        if($13!=".") {
            ensids[peak]=ensids[peak]?ensids[peak]";"$4:$4;  # ensembl id
        }
    } END{
        # Print gene info mapping
        for(peak in ensids) {
            print peak, ensids[peak]
        }
    }' > "$OUTPUT_DIR/temp_gene_info.txt"

# Run annotateBed for feature counting
echo "Running annotateBed for feature counting..."
bedtools annotate -counts \
    -i "$PEAKS" \
    -files $FEATURE_FILES \
    > "$OUTPUT_DIR/temp_feature_counts.bed"

# Combine feature counts with gene information
echo "Combining annotations..."
awk -v OFS='\t' '
    # Read gene info into arrays
    FILENAME == ARGV[1] {
        split($1, coords, "_");
        ens_info[coords[1],coords[2],coords[3]]=$2;
        next;
    }
    # Process feature counts and add gene info
    {
        ens_ids = ens_info[$1,$2,$3];
        if (!ens_ids) ens_ids = ".";
        print $0, ens_ids;
    }' "$OUTPUT_DIR/temp_gene_info.txt" "$OUTPUT_DIR/temp_feature_counts.bed" \
    > "$OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks.bed"

# Create R script for analysis
echo "Creating R script for analysis..."
cat > "$OUTPUT_DIR/analyze_features.R" << 'EOF'
# Read and analyze feature counts
library(tidyverse)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript analyze_features.R annotated_peaks.bed sample_id")
}

input_file <- args[1]
sample_id <- args[2]

# Read the feature counts file
data <- read.table(input_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Define feature names (matching the order in the shell script)
feature_names <- c("CDS", "Exons", "5'UTR", "Genes", "Introns", "3'UTR", "Transcripts")

# Rename columns
names(data) <- c("chr", "start", "end", "name", "score", "strand", 
                paste0(feature_names, "_count"), "ensembl_ids")

# Calculate summary statistics
summary_stats <- data %>%
  summarise(across(ends_with("_count"), 
                  list(
                    total = ~sum(. > 0),
                    percent = ~mean(. > 0) * 100
                  ))) %>%
  pivot_longer(everything(),
              names_to = c("feature", "stat"),
              names_pattern = "(.+)_count_(.+)",
              values_to = "value")

# Create summary table
summary_wide <- summary_stats %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    feature = str_remove(feature, "_count"),
    percent = round(percent, 2)
  )

# Add sample ID to the summary
summary_wide$sample_id <- sample_id

# Write results
write.csv(summary_wide, 
          file = file.path(dirname(input_file), 
                          paste0(sample_id, "_feature_overlap_summary.csv")),
          row.names = FALSE)

# Create visualization
pdf(file.path(dirname(input_file), 
              paste0(sample_id, "_feature_overlap_plot.pdf")), 
    width = 10, height = 6)

# Bar plot of overlap percentages
ggplot(summary_wide, aes(x = reorder(feature, -percent), y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", percent)), 
            vjust = -0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Feature Type",
       y = "Percentage of Peaks Overlapping",
       title = paste("Distribution of Peak Overlaps Across Genomic Features -", sample_id),
       subtitle = paste("Total peaks analyzed:", nrow(data)))

dev.off()

# Print summary
cat(sprintf("\nFeature Overlap Summary for %s:\n", sample_id))
print(summary_wide)

# Gene annotation summary
gene_summary <- data %>%
  summarise(
    total_peaks = n(),
    peaks_with_genes = sum(ensembl_ids != "."),
    percent_with_genes = mean(ensembl_ids != ".") * 100
  )

cat("\nGene Annotation Summary:\n")
print(gene_summary)
EOF

# Run R script
echo "Running R analysis..."
Rscript "$OUTPUT_DIR/analyze_features.R" \
    "$OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks.bed" \
    "$SAMPLE_ID"

# Clean up temporary files
rm "$OUTPUT_DIR/temp_gene_info.txt" "$OUTPUT_DIR/temp_feature_counts.bed"

# Print results
echo -e "\nAnnotation complete!"
echo "Results written to:"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks.bed"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_feature_overlap_summary.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_feature_overlap_plot.pdf"

# Display summary if available
if [ -f "$OUTPUT_DIR/${SAMPLE_ID}_feature_overlap_summary.csv" ]; then
    echo -e "\nFeature overlap summary:"
    cat "$OUTPUT_DIR/${SAMPLE_ID}_feature_overlap_summary.csv"
fi 