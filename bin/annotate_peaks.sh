#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -p PEAKS_BED -f FEATURES_DIR -o OUTPUT_DIR [-H]"
    echo "Annotate peaks with genomic features using bedtools annotate"
    echo ""
    echo "Arguments:"
    echo "  -p PEAKS_BED    Input peaks BED file (from Piranha)"
    echo "  -f FEATURES_DIR Directory containing feature BED files"
    echo "  -o OUTPUT_DIR   Output directory"
    echo "  -H             Only annotate host regions (chromosomes starting with 'chr')"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hp:f:o:H" opt; do
    case $opt in
        h) usage ;;
        p) PEAKS="$OPTARG" ;;
        f) FEATURES_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        H) HOST_ONLY=true ;;
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

# Filter peaks for host regions if requested
if [ "$HOST_ONLY" = true ]; then
    echo "Filtering for host regions (chr* chromosomes)..."
    TOTAL_PEAKS=$(wc -l < "$PEAKS")
    awk '$1 ~ /^chr/' "$PEAKS" > "$OUTPUT_DIR/host_peaks.bed"
    PEAKS="$OUTPUT_DIR/host_peaks.bed"
    HOST_PEAKS=$(wc -l < "$PEAKS")
    echo "Total peaks: $TOTAL_PEAKS"
    echo "Peaks in host regions: $HOST_PEAKS"
    echo "Peaks filtered out: $((TOTAL_PEAKS - HOST_PEAKS))"
    
    if [ "$HOST_PEAKS" -eq 0 ]; then
        echo "Error: No peaks found in host regions"
        rm "$OUTPUT_DIR/host_peaks.bed"
        exit 1
    fi
fi

# List of available feature types
FEATURES=("exons" "five_prime_utr" "genes" "introns" "three_prime_utr")

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

# Get gene annotations
echo "Getting gene annotations..."
bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/genes.bed" -wao | \
    awk 'BEGIN{OFS="\t"} {
        peak=$1"_"$2"_"$3;
        if($4!=".") {  # Assuming gene_id is in column 4
            genes[peak]=genes[peak]?genes[peak]";"$4:$4;  # gene_id
        }
    } END{
        for(peak in genes) {
            print peak, genes[peak];
        }
    }' > "$OUTPUT_DIR/temp_gene_info.txt"

# Get feature overlaps
echo "Getting feature overlaps..."
for feature in "${FEATURES[@]}"; do
    if [ "$feature" != "genes" ]; then
        bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/${feature}.bed" -wao | \
            awk -v feat="$feature" 'BEGIN{OFS="\t"} {
                peak=$1"_"$2"_"$3;
                if($4!=".") {  # Assuming feature_id is in column 4
                    features[peak]=features[peak]?features[peak]";"feat:feat;
                }
            } END{
                for(peak in features) {
                    print peak, features[peak];
                }
            }' >> "$OUTPUT_DIR/temp_features_info.txt"
    fi
done

# Run annotateBed for feature counting
echo "Running annotateBed for feature counting..."
bedtools annotate -counts \
    -i "$PEAKS" \
    -files $FEATURE_FILES \
    > "$OUTPUT_DIR/temp_feature_counts.bed"

# Combine all annotations
echo "Combining annotations..."
awk -v OFS='\t' '
    # Read gene info into arrays
    FILENAME == ARGV[1] {
        split($1, coords, "_");
        gene_ids[coords[1],coords[2],coords[3]]=$2;
        next;
    }
    # Read feature info into arrays
    FILENAME == ARGV[2] {
        split($1, coords, "_");
        feat_info[coords[1],coords[2],coords[3]]=$2;
        next;
    }
    # Process feature counts and add annotations
    {
        gene_id = gene_ids[$1,$2,$3];
        if (!gene_id) gene_id = ".";
        
        features = feat_info[$1,$2,$3];
        if (!features) features = ".";
        
        # Print original columns plus annotations
        print $0, gene_id, features;
    }' "$OUTPUT_DIR/temp_gene_info.txt" \
       "$OUTPUT_DIR/temp_features_info.txt" \
       "$OUTPUT_DIR/temp_feature_counts.bed" \
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
feature_names <- c("Exons", "5'UTR", "Genes", "Introns", "3'UTR")

# Rename columns
names(data) <- c("chr", "start", "end", "name", "score", "strand", 
                paste0(feature_names, "_count"), "gene_ids", "overlapping_features")

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
    peaks_with_genes = sum(gene_ids != "."),
    percent_with_genes = mean(gene_ids != ".") * 100
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
rm "$OUTPUT_DIR/temp_gene_info.txt" "$OUTPUT_DIR/temp_features_info.txt" "$OUTPUT_DIR/temp_feature_counts.bed"
if [ "$HOST_ONLY" = true ]; then
    rm "$OUTPUT_DIR/host_peaks.bed"
fi

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