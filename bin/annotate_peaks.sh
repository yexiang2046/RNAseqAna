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

# Get gene annotations with strand information
echo "Getting gene annotations..."
bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/genes.bed" -wao | \
    awk 'BEGIN{OFS="\t"; FS="\t"} {
        peak=$1"_"$2"_"$3;
        if($10!=".") {  # Use gene ID from column 10 (bedtools intersect output)
            genes[peak]=$10;  # Use gene ID
            strands[peak]=$6;  # Store strand information from the gene file
        }
    } END{
        for(peak in genes) {
            print peak, genes[peak], strands[peak];
        }
    }' > "$OUTPUT_DIR/temp_gene_info.txt"

# Get feature overlaps with priority (exons > UTRs > introns > genes > intergenic)
echo "Getting feature overlaps..."
# First process exons separately to ensure they take highest priority
echo "Processing exons..."
bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/exons.bed" -wao | \
    awk 'BEGIN{OFS="\t"} {
        peak=$1"_"$2"_"$3;
        overlap_len = $(NF);
        
        if(overlap_len > 0) {
            # Store overlap information with maximum overlap
            if(!(peak in max_overlap) || overlap_len > max_overlap[peak]) {
                max_overlap[peak] = overlap_len;
                features[peak] = "exons";
                # Extract gene ID from column 10 (bedtools intersect output)
                feature_ids[peak] = $10;
                # Store strand information from the feature file (column 12)
                feature_strands[peak] = $12;
            }
        }
    } END{
        for(peak in features) {
            print peak, feature_ids[peak], features[peak], max_overlap[peak], feature_strands[peak];
        }
    }' > "$OUTPUT_DIR/temp_exons_info.txt"

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
                    # Extract gene ID from column 10 (bedtools intersect output)
                    feature_ids[peak] = $10;
                    # Store strand information from the feature file (column 12)
                    feature_strands[peak] = $12;
                    seen[peak] = 1;
                }
            }
        } END{
            for(peak in features) {
                print peak, feature_ids[peak], features[peak], max_overlap[peak], feature_strands[peak];
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
                # Extract gene ID from column 10 (bedtools intersect output)
                feature_ids[peak] = $10;
                # Store strand information from the feature file (column 12)
                feature_strands[peak] = $12;
                seen[peak] = 1;
            }
        }
    } END{
        for(peak in features) {
            print peak, feature_ids[peak], features[peak], max_overlap[peak], feature_strands[peak];
        }
    }' > "$OUTPUT_DIR/temp_introns_info.txt"

# Combine feature files and add intergenic peaks
echo "Identifying intergenic peaks..."
cat "$OUTPUT_DIR"/temp_*_info.txt > "$OUTPUT_DIR/temp_all_features.txt"

# Find peaks without features and mark them as intergenic
bedtools intersect -a "$PEAKS" -b "$FEATURES_DIR/genes.bed" -v | \
    awk 'BEGIN{OFS="\t"} {
        peak=$1"_"$2"_"$3;
        print peak, ".", "intergenic", "0", ".";
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
                feature_strands[peak]=$5;  # Store strand information
                seen[peak]=1;
            }
        }
        END{
            for(peak in features) {
                print peak, feature_ids[peak], features[peak], feature_strands[peak];
            }
        }' | sort -k1,1 > "$OUTPUT_DIR/temp_features_info.txt"

# Combine all annotations
echo "Combining annotations..."
awk -v OFS='\t' '
    # Read feature info into arrays
    FILENAME == ARGV[1] {
        split($1, coords, "_");
        key = coords[1]"_"coords[2]"_"coords[3];
        feat_id[key]=$2;
        feat_type[key]=$3;
        feat_strand[key]=$4;  # Store feature strand
        next;
    }
    # Process peaks and add annotations
    {
        key = $1"_"$2"_"$3;
        feature_id = feat_id[key];
        feature_type = feat_type[key];
        feature_strand = feat_strand[key];
        if (!feature_type) {
            feature_id = ".";
            feature_type = "intergenic";
            feature_strand = ".";
        }
        
        # Print original columns but replace column 6 with feature strand
        # Use key as unique identifier to avoid duplicates
        if (!seen[key]) {
            print $1, $2, $3, $4, $5, feature_strand, feature_id, feature_type, feature_strand;
            seen[key] = 1;
        }
    }' "$OUTPUT_DIR/temp_features_info.txt" \
       "$PEAKS" \
    > "$OUTPUT_DIR/temp_annotated_peaks.bed"

# Sort and deduplicate the final output
echo "Sorting and deduplicating final output..."
sort -k1,1 -k2,2n -k3,3n "$OUTPUT_DIR/temp_annotated_peaks.bed" | \
    awk '!seen[$1"_"$2"_"$3]++' > "$OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks.bed"

# Add gene name and type annotations
echo "Adding gene name and type annotations..."
awk -v OFS='\t' '
    # First read genes.bed to create gene ID to name and type mapping
    FILENAME == ARGV[1] {
        if($4 != ".") {  # If gene ID exists
            gene_names[$4] = $5;  # Use column 5 as gene name
            gene_types[$4] = $7;  # Use column 7 as gene type
        }
        next;
    }
    # Process the annotated peaks file
    {
        # Get gene ID from column 7
        gene_id = $7;
        
        # Get gene name and type from mapping
        gene_name = gene_names[gene_id];
        if (!gene_name) gene_name = ".";
        
        gene_type = gene_types[gene_id];
        if (!gene_type) gene_type = ".";
        
        # Print all columns plus gene name and type
        print $0, gene_name, gene_type;
    }' "$FEATURES_DIR/genes.bed" \
       "$OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks.bed" \
    > "$OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks_with_names.bed"

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

# Rename columns (based on the actual columns in the annotated peaks file)
names(data) <- c("chr", "start", "end", "name", "score", "strand", 
                "gene_id", "feature_id", "feature_type", "feature_strand",
                "gene_name", "gene_type")

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

# Calculate feature type distribution
feature_stats <- data %>%
  group_by(feature_type) %>%
  summarise(
    count = n(),
    percentage = n() / nrow(data) * 100
  ) %>%
  arrange(desc(count))

# Write feature statistics
write.csv(feature_stats,
          file = file.path(dirname(input_file),
                          paste0(sample_id, "_feature_stats.csv")),
          row.names = FALSE)

# Create feature type distribution plot
pdf(file.path(dirname(input_file),
              paste0(sample_id, "_feature_distribution.pdf")),
    width = 10, height = 6)

ggplot(feature_stats, aes(x = reorder(feature_type, -count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, percentage)),
            vjust = -0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Feature Type",
       y = "Number of Peaks",
       title = paste("Distribution of Peak Features -", sample_id))

dev.off()

# Calculate strand distribution
strand_stats <- data %>%
  filter(feature_strand != ".") %>%
  group_by(feature_strand) %>%
  summarise(
    count = n(),
    percentage = n() / nrow(data) * 100
  ) %>%
  arrange(desc(count))

# Write strand statistics
write.csv(strand_stats,
          file = file.path(dirname(input_file),
                          paste0(sample_id, "_strand_stats.csv")),
          row.names = FALSE)

# Create strand distribution plot
pdf(file.path(dirname(input_file),
              paste0(sample_id, "_strand_distribution.pdf")),
    width = 10, height = 6)

ggplot(strand_stats, aes(x = feature_strand, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, percentage)),
            vjust = -0.5) +
  theme_minimal() +
  labs(x = "Strand",
       y = "Number of Peaks",
       title = paste("Distribution of Peak Strands -", sample_id))

dev.off()

# Calculate gene annotation statistics
gene_stats <- data %>%
  summarise(
    total_peaks = n(),
    peaks_with_genes = sum(gene_id != "."),
    peaks_with_names = sum(gene_name != "."),
    peaks_with_types = sum(gene_type != "."),
    peaks_with_strand = sum(feature_strand != "."),
    percent_with_genes = mean(gene_id != ".") * 100,
    percent_with_names = mean(gene_name != ".") * 100,
    percent_with_types = mean(gene_type != ".") * 100,
    percent_with_strand = mean(feature_strand != ".") * 100
  )

# Write gene annotation statistics
write.csv(gene_stats,
          file = file.path(dirname(input_file),
                          paste0(sample_id, "_gene_stats.csv")),
          row.names = FALSE)

# Calculate gene type distribution
gene_type_stats <- data %>%
  filter(gene_type != ".") %>%
  group_by(gene_type) %>%
  summarise(
    count = n(),
    percentage = n() / nrow(data) * 100
  ) %>%
  arrange(desc(count))

# Write gene type statistics
write.csv(gene_type_stats,
          file = file.path(dirname(input_file),
                          paste0(sample_id, "_gene_type_stats.csv")),
          row.names = FALSE)

# Create gene type distribution plot
pdf(file.path(dirname(input_file),
              paste0(sample_id, "_gene_type_distribution.pdf")),
    width = 10, height = 6)

ggplot(gene_type_stats, aes(x = reorder(gene_type, -count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, percentage)),
            vjust = -0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene Type",
       y = "Number of Peaks",
       title = paste("Distribution of Peak Gene Types -", sample_id))

dev.off()

# Create top genes table
top_genes <- data %>%
  filter(gene_name != ".") %>%
  group_by(gene_name, gene_type, feature_strand) %>%
  summarise(
    peak_count = n(),
    feature_types = paste(unique(feature_type), collapse=", ")
  ) %>%
  arrange(desc(peak_count)) %>%
  head(20)

# Write top genes table
write.csv(top_genes,
          file = file.path(dirname(input_file),
                          paste0(sample_id, "_top_genes.csv")),
          row.names = FALSE)

# Print summary
cat(sprintf("\nFeature Distribution for %s:\n", sample_id))
print(feature_stats)

cat("\nPeak Length Statistics:\n")
print(length_stats)

cat("\nStrand Distribution:\n")
print(strand_stats)

cat("\nGene Annotation Summary:\n")
print(gene_stats)

cat("\nGene Type Distribution:\n")
print(gene_type_stats)

cat("\nTop 20 Genes with Most Peaks:\n")
print(top_genes)
EOF

# Run R script
echo "Running R analysis..."
Rscript "$OUTPUT_DIR/analyze_features.R" \
    "$OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks_with_names.bed" \
    "$SAMPLE_ID"

# Clean up temporary files
rm "$OUTPUT_DIR/temp_gene_info.txt" "$OUTPUT_DIR/temp_features_info.txt" "$OUTPUT_DIR/temp_all_features.txt" "$OUTPUT_DIR/temp_intergenic.txt" "$OUTPUT_DIR/temp_annotated_peaks.bed"
if [ "$HOST_ONLY" = true ]; then
    rm "$OUTPUT_DIR/host_peaks.bed"
fi

# Print results
echo -e "\nAnnotation complete!"
echo "Results written to:"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_annotated_peaks_with_names.bed"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_feature_stats.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_feature_distribution.pdf"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_gene_stats.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_gene_type_stats.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_gene_type_distribution.pdf"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_top_genes.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_peak_length_stats.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_peak_length_distribution.pdf"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_strand_stats.csv"
echo "- $OUTPUT_DIR/${SAMPLE_ID}_strand_distribution.pdf"

# Display summary if available
if [ -f "$OUTPUT_DIR/${SAMPLE_ID}_feature_stats.csv" ]; then
    echo -e "\nFeature statistics:"
    cat "$OUTPUT_DIR/${SAMPLE_ID}_feature_stats.csv"
fi

if [ -f "$OUTPUT_DIR/${SAMPLE_ID}_peak_length_stats.csv" ]; then
    echo -e "\nPeak length statistics:"
    cat "$OUTPUT_DIR/${SAMPLE_ID}_peak_length_stats.csv"
fi

if [ -f "$OUTPUT_DIR/${SAMPLE_ID}_strand_stats.csv" ]; then
    echo -e "\nStrand statistics:"
    cat "$OUTPUT_DIR/${SAMPLE_ID}_strand_stats.csv"
fi 