#!/bin/bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i BED_DIR -o OUTPUT_DIR [-m MIN_OVERLAP]"
    echo "Find overlapping regions between BED files"
    echo ""
    echo "Arguments:"
    echo "  -i BED_DIR      Directory containing BED files"
    echo "  -o OUTPUT_DIR   Output directory"
    echo "  -m MIN_OVERLAP  Minimum overlap required (default: 1bp)"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hi:o:m:" opt; do
    case $opt in
        h) usage ;;
        i) BED_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) MIN_OVERLAP="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$BED_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Set default minimum overlap if not specified
MIN_OVERLAP=${MIN_OVERLAP:-1}

# Check if input directory exists
if [ ! -d "$BED_DIR" ]; then
    echo "Error: BED directory does not exist: $BED_DIR"
    exit 1
fi

# Check if bedtools is installed
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR/pairwise" "$OUTPUT_DIR/summary"

# Get list of BED files
BED_FILES=($(ls "$BED_DIR"/*.bed 2>/dev/null))
NUM_FILES=${#BED_FILES[@]}

if [ $NUM_FILES -lt 2 ]; then
    echo "Error: At least 2 BED files are required. Found: $NUM_FILES"
    exit 1
fi

echo "Found $NUM_FILES BED files"

# Create summary file
SUMMARY_FILE="$OUTPUT_DIR/summary/overlap_summary.txt"
echo "Pairwise Overlap Analysis Summary" > "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo "Minimum overlap required: ${MIN_OVERLAP}bp" >> "$SUMMARY_FILE"
echo "----------------------------------------" >> "$SUMMARY_FILE"

# Create CSV summary
CSV_SUMMARY="$OUTPUT_DIR/summary/overlap_summary.csv"
echo "File1,File2,Total_Regions_1,Total_Regions_2,Overlapping_Regions,Overlap_Percentage" > "$CSV_SUMMARY"

# Function to get base filename without extension
get_basename() {
    basename "$1" .bed
}

# Process each pair of files
for ((i=0; i<$NUM_FILES; i++)); do
    for ((j=i+1; j<$NUM_FILES; j++)); do
        file1="${BED_FILES[$i]}"
        file2="${BED_FILES[$j]}"
        
        base1=$(get_basename "$file1")
        base2=$(get_basename "$file2")
        
        echo "Processing: $base1 vs $base2"
        
        # Count total regions in each file
        total1=$(wc -l < "$file1")
        total2=$(wc -l < "$file2")
        
        # Find overlaps
        overlap_file="$OUTPUT_DIR/pairwise/${base1}_vs_${base2}_overlap.bed"
        bedtools intersect -a "$file1" -b "$file2" -wa -wb \
            -f ${MIN_OVERLAP} -r > "$overlap_file"
        
        # Count unique overlapping regions
        overlaps=$(cut -f1-3 "$overlap_file" | sort -u | wc -l)
        
        # Calculate percentage (relative to file with fewer regions)
        min_regions=$(($total1 < $total2 ? $total1 : $total2))
        if [ $min_regions -eq 0 ]; then
            percentage=0
        else
            percentage=$(echo "scale=2; ($overlaps * 100) / $min_regions" | bc)
        fi
        
        # Add to summary
        echo -e "\n$base1 vs $base2:" >> "$SUMMARY_FILE"
        echo "Total regions in $base1: $total1" >> "$SUMMARY_FILE"
        echo "Total regions in $base2: $total2" >> "$SUMMARY_FILE"
        echo "Overlapping regions: $overlaps" >> "$SUMMARY_FILE"
        echo "Overlap percentage: ${percentage}%" >> "$SUMMARY_FILE"
        
        # Add to CSV
        echo "$base1,$base2,$total1,$total2,$overlaps,$percentage" >> "$CSV_SUMMARY"
        
        # Create BED file with overlap statistics
        overlap_stats="$OUTPUT_DIR/pairwise/${base1}_vs_${base2}_stats.bed"
        awk -v OFS='\t' '
        {
            # Combine coordinates to get overlap region
            start = $2 > $8 ? $2 : $8
            end = $3 < $9 ? $3 : $9
            overlap_size = end - start
            
            # Print overlap region with statistics
            print $1, start, end, \
                  "overlap_size=" overlap_size, \
                  "file1_region=" $1 ":" $2 "-" $3, \
                  "file2_region=" $7 ":" $8 "-" $9
        }' "$overlap_file" > "$overlap_stats"
    done
done

# Create R script for visualization
cat > "$OUTPUT_DIR/summary/plot_overlaps.R" << 'EOF'
library(ggplot2)
library(reshape2)

# Read the data
data <- read.csv("overlap_summary.csv")

# Create heatmap matrix
matrix_data <- reshape2::acast(data, File1~File2, value.var="Overlap_Percentage")

# Create PDF
pdf("overlap_heatmap.pdf", width=10, height=8)

# Plot heatmap
heatmap(matrix_data, 
        main="Overlap Percentage Heatmap",
        col=colorRampPalette(c("white", "steelblue"))(100),
        margins=c(10,10))

# Create barplot of overlap percentages
ggplot(data, aes(x=paste(File1, "vs", File2), y=Overlap_Percentage)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x="File Pairs", y="Overlap Percentage (%)",
         title="Pairwise Overlap Percentages")

dev.off()
EOF

# Run R script if R is available
if command -v Rscript &> /dev/null; then
    echo "Creating visualizations..."
    cd "$OUTPUT_DIR/summary" && Rscript plot_overlaps.R
fi

echo "Analysis complete!"
echo "Results written to: $OUTPUT_DIR"
echo "See summary file: $SUMMARY_FILE"
echo "See CSV summary: $CSV_SUMMARY"
if [ -f "$OUTPUT_DIR/summary/overlap_heatmap.pdf" ]; then
    echo "See visualizations: $OUTPUT_DIR/summary/overlap_heatmap.pdf"
fi 