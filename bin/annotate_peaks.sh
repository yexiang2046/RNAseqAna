#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -p PEAKS_BED -f FEATURES_DIR -o OUTPUT_DIR [-H] [-m MIN_LENGTH] [-M MAX_LENGTH]"
    echo "Annotate peaks with genomic features using bedtools annotate"
    echo ""
    echo "Arguments:"
    echo "  -p PEAKS_BED    Input peaks BED file (from Piranha)"
    echo "  -f FEATURES_DIR Directory containing feature BED files"
    echo "  -o OUTPUT_DIR   Output directory"
    echo "  -H             Only annotate host regions (chromosomes starting with 'chr')"
    echo "  -m MIN_LENGTH  Minimum peak length to include (default: no minimum)"
    echo "  -M MAX_LENGTH  Maximum peak length to include (default: no maximum)"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hp:f:o:Hm:M:" opt; do
    case $opt in
        h) usage ;;
        p) PEAKS="$OPTARG" ;;
        f) FEATURES_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        H) HOST_ONLY=true ;;
        m) MIN_LENGTH="$OPTARG" ;;
        M) MAX_LENGTH="$OPTARG" ;;
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

# Create temporary working file
cp "$PEAKS" "$OUTPUT_DIR/temp_peaks.bed"
PEAKS="$OUTPUT_DIR/temp_peaks.bed"

# Filter peaks by length if specified
if [ ! -z "$MIN_LENGTH" ] || [ ! -z "$MAX_LENGTH" ]; then
    echo "Filtering peaks by length..."
    TOTAL_PEAKS=$(wc -l < "$PEAKS")
    
    # Construct awk condition for length filtering
    LENGTH_CONDITION=""
    if [ ! -z "$MIN_LENGTH" ]; then
        LENGTH_CONDITION="(\$3-\$2) >= $MIN_LENGTH"
    fi
    if [ ! -z "$MAX_LENGTH" ]; then
        if [ ! -z "$LENGTH_CONDITION" ]; then
            LENGTH_CONDITION="$LENGTH_CONDITION && "
        fi
        LENGTH_CONDITION="${LENGTH_CONDITION}(\$3-\$2) <= $MAX_LENGTH"
    fi
    
    # Apply length filter
    awk -v OFS="\t" "$LENGTH_CONDITION" "$PEAKS" > "$OUTPUT_DIR/length_filtered_peaks.bed"
    mv "$OUTPUT_DIR/length_filtered_peaks.bed" "$PEAKS"
    
    FILTERED_PEAKS=$(wc -l < "$PEAKS")
    echo "Total peaks: $TOTAL_PEAKS"
    echo "Peaks after length filtering: $FILTERED_PEAKS"
    echo "Peaks filtered out: $((TOTAL_PEAKS - FILTERED_PEAKS))"
    
    if [ "$FILTERED_PEAKS" -eq 0 ]; then
        echo "Error: No peaks remain after length filtering"
        rm "$PEAKS"
        exit 1
    fi
fi

# Filter peaks for host regions if requested
if [ "$HOST_ONLY" = true ]; then
    echo "Filtering for host regions (chr* chromosomes)..."
    TOTAL_PEAKS=$(wc -l < "$PEAKS")
    awk '$1 ~ /^chr/' "$PEAKS" > "$OUTPUT_DIR/host_peaks.bed"
    mv "$OUTPUT_DIR/host_peaks.bed" "$PEAKS"
    HOST_PEAKS=$(wc -l < "$PEAKS")
    echo "Total peaks: $TOTAL_PEAKS"
    echo "Peaks in host regions: $HOST_PEAKS"
    echo "Peaks filtered out: $((TOTAL_PEAKS - HOST_PEAKS))"
    
    if [ "$HOST_PEAKS" -eq 0 ]; then
        echo "Error: No peaks found in host regions"
        rm "$PEAKS"
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
    awk 'BEGIN{OFS="\t"; FS="\t"} {
        peak=$1"_"$2"_"$3;
        if($10!=".") {  # Use gene ID from column 10 (bedtools intersect output)
            genes[peak]=$10;  # Use gene ID
        }
    } END{
        for(peak in genes) {
            print peak, genes[peak];
        }
    }' > "$OUTPUT_DIR/temp_gene_info.txt"

# Get feature overlaps with priority (exons > UTRs > introns > genes > intergenic)
echo "Getting feature overlaps..."
# First process exons separately to ensure they take highest priority
echo "Processing exons..."
echo "Checking input files..."
echo "Number of peaks in input file:"
wc -l "$PEAKS"
echo "Number of exons in feature file:"
wc -l "$FEATURES_DIR/exons.bed"

# Clean input files to ensure proper tab formatting
echo "Cleaning input files..."
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6}' "$PEAKS" > "$OUTPUT_DIR/clean_peaks.bed"
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6}' "$FEATURES_DIR/exons.bed" > "$OUTPUT_DIR/clean_exons.bed"

echo "First few lines of cleaned peaks file:"
head -n 3 "$OUTPUT_DIR/clean_peaks.bed"
echo "First few lines of cleaned exons file:"
head -n 3 "$OUTPUT_DIR/clean_exons.bed"

echo "Running bedtools intersect for exons..."
bedtools intersect -a "$OUTPUT_DIR/clean_peaks.bed" -b "$OUTPUT_DIR/clean_exons.bed" -wao > "$OUTPUT_DIR/debug_exon_overlaps.txt"

echo "Checking if bedtools intersect produced output..."
if [ ! -s "$OUTPUT_DIR/debug_exon_overlaps.txt" ]; then
    echo "Warning: bedtools intersect produced no output!"
    echo "Checking if bedtools intersect command works with a test case..."
    # Create a small test file
    echo -e "chr1\t100\t200\tpeak1\t0\t+" > "$OUTPUT_DIR/test_peaks.bed"
    echo -e "chr1\t150\t250\texon1\t0\t+" > "$OUTPUT_DIR/test_exons.bed"
    bedtools intersect -a "$OUTPUT_DIR/test_peaks.bed" -b "$OUTPUT_DIR/test_exons.bed" -wao
    rm "$OUTPUT_DIR/test_peaks.bed" "$OUTPUT_DIR/test_exons.bed"
fi

echo "Checking exon overlaps..."
awk 'BEGIN{OFS="\t"; FS="\t"} {
    # For bedtools intersect -wao output:
    # $1-$3: peak coordinates
    # $4-$6: exon coordinates
    # $7: gene ID from exon file (column 4 of original exon.bed)
    # $8: score
    # $9: strand
    # $NF: overlap length
    peak=$1"_"$2"_"$3;
    overlap_len = $(NF);
    
    print "Processing peak:", peak, "overlap length:", overlap_len > "/dev/stderr";
    
    if(overlap_len > 0) {
        # Store overlap information with maximum overlap
        if(!(peak in max_overlap) || overlap_len > max_overlap[peak]) {
            max_overlap[peak] = overlap_len;
            features[peak] = "exons";
            # Extract gene ID from column 7 (bedtools intersect output)
            feature_ids[peak] = $7;
            print "Found exon overlap for peak:", peak, "gene ID:", $7 > "/dev/stderr";
        }
    }
} END{
    print "Total peaks with exon overlaps:", length(features) > "/dev/stderr";
    for(peak in features) {
        print peak, feature_ids[peak], features[peak], max_overlap[peak];
    }
}' "$OUTPUT_DIR/debug_exon_overlaps.txt" > "$OUTPUT_DIR/temp_exons_info.txt"

# Clean up temporary files
rm "$OUTPUT_DIR/clean_peaks.bed" "$OUTPUT_DIR/clean_exons.bed"

# Then process UTRs
for feature in "five_prime_utr" "three_prime_utr"; do
    echo "Processing $feature..."
    bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/${feature}.bed" -wao | \
        awk -v feat="$feature" 'BEGIN{OFS="\t"} {
            peak=$1"_"$2"_"$3;
            overlap_len = $(NF);
            
            if(overlap_len > 0) {
                # Only store if no exon overlap exists
                if(!(peak in seen)) {
                    max_overlap[peak] = overlap_len;
                    features[peak] = feat;
                    # Extract gene ID from column 7 (bedtools intersect output)
                    feature_ids[peak] = $7;
                    seen[peak] = 1;
                }
            }
        } END{
            for(peak in features) {
                print peak, feature_ids[peak], features[peak], max_overlap[peak];
            }
        }' > "$OUTPUT_DIR/temp_${feature}_info.txt"
done

# Then process introns
echo "Processing introns..."
bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/introns.bed" -wao | \
    awk 'BEGIN{OFS="\t"} {
        peak=$1"_"$2"_"$3;
        overlap_len = $(NF);
        
        if(overlap_len > 0) {
            # Only store if no exon or UTR overlap exists
            if(!(peak in seen)) {
                max_overlap[peak] = overlap_len;
                features[peak] = "introns";
                # Extract gene ID from column 7 (bedtools intersect output)
                feature_ids[peak] = $7;
                seen[peak] = 1;
            }
        }
    } END{
        for(peak in features) {
            print peak, feature_ids[peak], features[peak], max_overlap[peak];
        }
    }' > "$OUTPUT_DIR/temp_introns_info.txt"

# Combine feature files and add intergenic peaks
echo "Identifying intergenic peaks..."
cat "$OUTPUT_DIR"/temp_*_info.txt > "$OUTPUT_DIR/temp_all_features.txt"

# Find peaks without features and mark them as intergenic
bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/genes.bed" -v | \
    awk 'BEGIN{OFS="\t"} {
        peak=$1"_"$2"_"$3;
        print peak, ".", "intergenic", "0";
    }' > "$OUTPUT_DIR/temp_intergenic.txt"

# Sort and prioritize features (exons > UTRs > introns > genes > intergenic)
cat "$OUTPUT_DIR/temp_all_features.txt" "$OUTPUT_DIR/temp_intergenic.txt" | \
    awk 'BEGIN{OFS="\t"; 
             # Define feature priority
             priority["exons"]=1;
             priority["five_prime_utr"]=2;
             priority["three_prime_utr"]=2;
             priority["introns"]=3;
             priority["genes"]=4;
             priority["intergenic"]=5;
        }
        {
            peak=$1;
            feat_type=$3;
            curr_priority = priority[feat_type];
            if(curr_priority == "") { curr_priority = 999; }  # Handle unknown features
            
            if(!(peak in seen) || curr_priority < min_priority[peak]) {
                features[peak]=feat_type;
                feature_ids[peak]=$2;
                min_priority[peak] = curr_priority;
                seen[peak]=1;
            }
        }
        END{
            for(peak in features) {
                print peak, feature_ids[peak], features[peak];
            }
        }' | sort -k1,1 > "$OUTPUT_DIR/temp_features_info.txt"

# Combine all annotations
echo "Combining annotations..."
awk -v OFS='\t' '
    # Read gene info into arrays
    FILENAME == ARGV[1] {
        split($1, coords, "_");
        key = coords[1]"_"coords[2]"_"coords[3];
        gene_ids[key]=$2;
        next;
    }
    # Read feature info into arrays
    FILENAME == ARGV[2] {
        split($1, coords, "_");
        key = coords[1]"_"coords[2]"_"coords[3];
        feat_id[key]=$2;
        feat_type[key]=$3;
        next;
    }
    # Process feature counts and add annotations
    {
        key = $1"_"$2"_"$3;
        gene_id = gene_ids[key];
        if (!gene_id) gene_id = ".";
        
        feature_id = feat_id[key];
        feature_type = feat_type[key];
        if (!feature_type) {
            feature_id = ".";
            feature_type = "intergenic";
        }
        
        # Print original columns plus annotations
        print $0, gene_id, feature_id, feature_type;
    }' "$OUTPUT_DIR/temp_gene_info.txt" \
       "$OUTPUT_DIR/temp_features_info.txt" \
       "$PEAKS" \
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

# Calculate peak lengths and statistics
data <- data %>%
  mutate(peak_length = end - start)

length_stats <- data %>%
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
          file = file.path(dirname(input_file),
                          paste0(sample_id, "_peak_length_stats.csv")),
          row.names = FALSE)

# Create length distribution plot
pdf(file.path(dirname(input_file),
              paste0(sample_id, "_peak_length_distribution.pdf")),
    width = 10, height = 6)

ggplot(data, aes(x = peak_length)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(x = "Peak Length (bp)",
       y = "Count",
       title = paste("Distribution of Peak Lengths -", sample_id),
       subtitle = sprintf("Mean: %.1f bp, Median: %.1f bp",
                        length_stats$mean_length,
                        length_stats$median_length))

dev.off()

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

# Print length statistics
cat("\nPeak Length Statistics:\n")
print(length_stats)

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
rm "$OUTPUT_DIR/temp_gene_info.txt" "$OUTPUT_DIR/temp_features_info.txt" "$OUTPUT_DIR/temp_all_features.txt" "$OUTPUT_DIR/temp_intergenic.txt"
if [ "$HOST_ONLY" = true ]; then
    rm "$OUTPUT_DIR/host_peaks.bed"
fi

# Print results
echo -e "\nAnnotation complete!"
echo "Results written to:"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks.bed"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_feature_overlap_summary.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_feature_overlap_plot.pdf"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_peak_length_stats.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_peak_length_distribution.pdf"

# Display summary if available
if [ -f "$OUTPUT_DIR/${SAMPLE_ID}_feature_overlap_summary.csv" ]; then
    echo -e "\nFeature overlap summary:"
    cat "$OUTPUT_DIR/${SAMPLE_ID}_feature_overlap_summary.csv"
fi

if [ -f "$OUTPUT_DIR/${SAMPLE_ID}_peak_length_stats.csv" ]; then
    echo -e "\nPeak length statistics:"
    cat "$OUTPUT_DIR/${SAMPLE_ID}_peak_length_stats.csv"
fi 