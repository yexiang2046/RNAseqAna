#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i BED_DIR -o OUTPUT_DIR [-m MIN_OVERLAP]"
    echo "Find overlapping regions between BED files (pairwise and three-way)"
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
mkdir -p "$OUTPUT_DIR/pairwise" "$OUTPUT_DIR/three_way" "$OUTPUT_DIR/summary"

# Get list of BED files
BED_FILES=($(ls "$BED_DIR"/*.bed 2>/dev/null))
NUM_FILES=${#BED_FILES[@]}

if [ $NUM_FILES -lt 2 ]; then
    echo "Error: At least 2 BED files are required. Found: $NUM_FILES"
    exit 1
fi

echo "Found $NUM_FILES BED files"

# Create summary files
SUMMARY_FILE="$OUTPUT_DIR/summary/overlap_summary.txt"
echo "Overlap Analysis Summary" > "$SUMMARY_FILE"
echo "Date: $(date)" >> "$SUMMARY_FILE"
echo "Minimum overlap required: ${MIN_OVERLAP}bp" >> "$SUMMARY_FILE"
echo "----------------------------------------" >> "$SUMMARY_FILE"

# Create CSV summaries
PAIRWISE_CSV="$OUTPUT_DIR/summary/pairwise_overlap_summary.csv"
THREE_WAY_CSV="$OUTPUT_DIR/summary/three_way_overlap_summary.csv"
echo "File1,File2,Total_Regions_1,Total_Regions_2,Overlapping_Regions,Overlap_Percentage" > "$PAIRWISE_CSV"
echo "File1,File2,File3,Total_Regions_1,Total_Regions_2,Total_Regions_3,Common_Regions,Overlap_Percentage" > "$THREE_WAY_CSV"

# Function to get base filename without extension
get_basename() {
    basename "$1" .bed
}

echo "Processing pairwise overlaps..."
# Process each pair of files
for ((i=0; i<$NUM_FILES; i++)); do
    for ((j=i+1; j<$NUM_FILES; j++)); do
        file1="${BED_FILES[$i]}"
        file2="${BED_FILES[$j]}"
        
        base1=$(get_basename "$file1")
        base2=$(get_basename "$file2")
        
        echo "Processing pair: $base1 vs $base2"
        
        # Count total regions in each file
        total1=$(wc -l < "$file1")
        total2=$(wc -l < "$file2")
        
        # Find overlaps and remove duplicates
        overlap_file="$OUTPUT_DIR/pairwise/${base1}_vs_${base2}_overlap.bed"
        bedtools intersect -s -a "$file1" -b "$file2" -wa -wb \
            -f ${MIN_OVERLAP} -r | \
            sort -k1,1 -k2,2n -k3,3n | \
            uniq > "$overlap_file"
        
        # Count unique overlapping regions
        overlaps=$(cut -f1-3 "$overlap_file" | sort -k1,1 -k2,2n -k3,3n | uniq | wc -l)
        
        # Calculate percentage (relative to file with fewer regions)
        min_regions=$(($total1 < $total2 ? $total1 : $total2))
        if [ $min_regions -eq 0 ]; then
            percentage=0
        else
            percentage=$(echo "scale=2; ($overlaps * 100) / $min_regions" | bc)
        fi
        
        # Add to summary
        echo -e "\nPairwise: $base1 vs $base2:" >> "$SUMMARY_FILE"
        echo "Total regions in $base1: $total1" >> "$SUMMARY_FILE"
        echo "Total regions in $base2: $total2" >> "$SUMMARY_FILE"
        echo "Overlapping regions: $overlaps" >> "$SUMMARY_FILE"
        echo "Overlap percentage: ${percentage}%" >> "$SUMMARY_FILE"
        
        # Add to CSV
        echo "$base1,$base2,$total1,$total2,$overlaps,$percentage" >> "$PAIRWISE_CSV"
        
        # Create BED file with overlap statistics
        overlap_stats="$OUTPUT_DIR/pairwise/${base1}_vs_${base2}_stats.bed"
        awk -v OFS='\t' '
        {
            start = $2 > $8 ? $2 : $8
            end = $3 < $9 ? $3 : $9
            overlap_size = end - start
            print $1, start, end, \
                  "overlap_size=" overlap_size, \
                  "file1_region=" $1 ":" $2 "-" $3, \
                  "file2_region=" $7 ":" $8 "-" $9
        }' "$overlap_file" > "$overlap_stats"
    done
done

echo "Processing three-way overlaps..."
# Process three-way overlaps
if [ $NUM_FILES -ge 3 ]; then
    for ((i=0; i<$NUM_FILES-2; i++)); do
        for ((j=i+1; j<$NUM_FILES-1; j++)); do
            for ((k=j+1; k<$NUM_FILES; k++)); do
                file1="${BED_FILES[$i]}"
                file2="${BED_FILES[$j]}"
                file3="${BED_FILES[$k]}"
                
                base1=$(get_basename "$file1")
                base2=$(get_basename "$file2")
                base3=$(get_basename "$file3")
                
                echo "Processing trio: $base1 vs $base2 vs $base3"
                
                # Count total regions
                total1=$(wc -l < "$file1")
                total2=$(wc -l < "$file2")
                total3=$(wc -l < "$file3")
                
                # Find three-way overlaps and remove duplicates
                overlap_file="$OUTPUT_DIR/three_way/${base1}_${base2}_${base3}_overlap.bed"
                temp_overlap="$OUTPUT_DIR/three_way/temp_overlap.bed"
                
                # First find overlap between first two files
                bedtools intersect -s -a "$file1" -b "$file2" -wa -wb \
                    -f ${MIN_OVERLAP} -r | \
                    sort -k1,1 -k2,2n -k3,3n | \
                    uniq > "$temp_overlap"
                
                # Then find overlap with third file and remove duplicates
                bedtools intersect -s -a "$temp_overlap" -b "$file3" -wa -wb \
                    -f ${MIN_OVERLAP} -r | \
                    sort -k1,1 -k2,2n -k3,3n | \
                    uniq > "$overlap_file"
                
                rm -f "$temp_overlap"
                
                # Count unique three-way overlapping regions
                common_regions=$(cut -f1-3 "$overlap_file" | sort -k1,1 -k2,2n -k3,3n | uniq | wc -l)
                
                # Calculate percentage (relative to file with fewest regions)
                min_regions=$(($total1 < $total2 ? $total1 : $total2))
                min_regions=$(($min_regions < $total3 ? $min_regions : $total3))
                if [ $min_regions -eq 0 ]; then
                    percentage=0
                else
                    percentage=$(echo "scale=2; ($common_regions * 100) / $min_regions" | bc)
                fi
                
                # Add to summary
                echo -e "\nThree-way: $base1 vs $base2 vs $base3:" >> "$SUMMARY_FILE"
                echo "Total regions in $base1: $total1" >> "$SUMMARY_FILE"
                echo "Total regions in $base2: $total2" >> "$SUMMARY_FILE"
                echo "Total regions in $base3: $total3" >> "$SUMMARY_FILE"
                echo "Common regions: $common_regions" >> "$SUMMARY_FILE"
                echo "Overlap percentage: ${percentage}%" >> "$SUMMARY_FILE"
                
                # Add to CSV
                echo "$base1,$base2,$base3,$total1,$total2,$total3,$common_regions,$percentage" >> "$THREE_WAY_CSV"
                
                # Create BED file with overlap statistics
                overlap_stats="$OUTPUT_DIR/three_way/${base1}_${base2}_${base3}_stats.bed"
                awk -v OFS='\t' '
                {
                    # Get the overlapping region coordinates
                    start = $2
                    end = $3
                    if ($8 > start) start = $8
                    if ($9 < end) end = $9
                    if ($14 > start) start = $14
                    if ($15 < end) end = $15
                    
                    overlap_size = end - start
                    
                    print $1, start, end, \
                          "overlap_size=" overlap_size, \
                          "file1_region=" $1 ":" $2 "-" $3, \
                          "file2_region=" $7 ":" $8 "-" $9, \
                          "file3_region=" $13 ":" $14 "-" $15
                }' "$overlap_file" > "$overlap_stats"
            done
        done
    done
fi

# Create R script for visualization
cat > "$OUTPUT_DIR/summary/plot_overlaps.R" << 'EOF'
library(ggplot2)
library(VennDiagram)
library(gridExtra)

# Function to create a color-transparent version of a color
transparent_color <- function(color, alpha=0.5) {
    rgb_col <- col2rgb(color)
    rgb(rgb_col[1], rgb_col[2], rgb_col[3], alpha=alpha*255, maxColorValue=255)
}

# Read the data
pairwise_data <- read.csv("pairwise_overlap_summary.csv")
three_way_data <- NULL
if (file.exists("three_way_overlap_summary.csv")) {
    three_way_data <- read.csv("three_way_overlap_summary.csv")
}

# Create PDF
pdf("overlap_visualization.pdf", width=12, height=10)

# Function to get non-overlapping counts for Venn diagram
get_venn_counts <- function(data) {
    unique_files <- unique(c(as.character(data$File1), as.character(data$File2)))
    total_counts <- sapply(unique_files, function(f) {
        rows <- data[data$File1 == f | data$File2 == f, ]
        if(length(rows) > 0) {
            return(rows[1, grep("Total_Regions", colnames(data))][1])
        }
        return(0)
    })
    names(total_counts) <- unique_files
    
    overlap_counts <- setNames(data$Overlapping_Regions, 
                             paste(data$File1, data$File2, sep="|"))
    
    return(list(totals=total_counts, overlaps=overlap_counts))
}

# Create Venn diagrams based on number of files
n_files <- length(unique(c(pairwise_data$File1, pairwise_data$File2)))

# Set up the plotting area
if (n_files >= 3 && !is.null(three_way_data)) {
    layout(matrix(c(1,2,3,3), 2, 2))
} else {
    layout(matrix(c(1,2), 1, 2))
}

# Colors for the Venn diagrams
venn_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")

# Create pairwise Venn diagrams
unique_pairs <- unique(paste(pairwise_data$File1, pairwise_data$File2))
par(mar=c(2,2,4,2))

# Split pairs into groups of 2 for visualization
n_pairs <- length(unique_pairs)
pairs_per_plot <- min(3, n_pairs)
pair_groups <- split(unique_pairs, ceiling(seq_along(unique_pairs)/pairs_per_plot))

for (pair_group in pair_groups) {
    # Create a new plot for each group of pairs
    plot.new()
    title("Pairwise Overlaps")
    
    for (i in seq_along(pair_group)) {
        pair <- strsplit(pair_group[i], " ")[[1]]
        pair_data <- pairwise_data[pairwise_data$File1 == pair[1] & 
                                  pairwise_data$File2 == pair[2], ]
        
        # Calculate positions for multiple Venn diagrams
        x_pos <- 0.2 + (i-1) * 0.3
        y_pos <- 0.5
        
        # Draw Venn diagram
        draw.pairwise.venn(
            area1 = pair_data$Total_Regions_1,
            area2 = pair_data$Total_Regions_2,
            cross.area = pair_data$Overlapping_Regions,
            category = c(pair[1], pair[2]),
            col = venn_colors[c(i*2-1, i*2)],
            fill = transparent_color(venn_colors[c(i*2-1, i*2)]),
            scaled = TRUE,
            euler.d = TRUE,
            sep.dist = 0.03,
            rotation.degree = 0,
            lwd = 2,
            cex = 0.8,
            cat.cex = 0.8,
            cat.dist = 0.05,
            ext.pos = 0,
            ext.dist = -0.05,
            ext.length = 0.85,
            ext.line.lwd = 2,
            ext.line.lty = "dashed"
        )
    }
}

# Create three-way Venn diagram if applicable
if (n_files >= 3 && !is.null(three_way_data)) {
    # Get a representative trio
    trio <- three_way_data[1,]
    
    plot.new()
    title("Three-way Overlaps")
    
    draw.triple.venn(
        area1 = trio$Total_Regions_1,
        area2 = trio$Total_Regions_2,
        area3 = trio$Total_Regions_3,
        n12 = sum(pairwise_data$Overlapping_Regions[
            pairwise_data$File1 == trio$File1 & 
            pairwise_data$File2 == trio$File2]),
        n23 = sum(pairwise_data$Overlapping_Regions[
            pairwise_data$File1 == trio$File2 & 
            pairwise_data$File2 == trio$File3]),
        n13 = sum(pairwise_data$Overlapping_Regions[
            pairwise_data$File1 == trio$File1 & 
            pairwise_data$File2 == trio$File3]),
        n123 = trio$Common_Regions,
        category = c(trio$File1, trio$File2, trio$File3),
        col = venn_colors[1:3],
        fill = sapply(venn_colors[1:3], transparent_color),
        scaled = TRUE,
        euler.d = TRUE,
        lwd = 2,
        cex = 0.8,
        cat.cex = 0.8,
        cat.col = venn_colors[1:3]
    )
}

# Create overlap percentage barplots
par(mar=c(10,4,4,2))
plot.new()

# Pairwise overlap barplot
p1 <- ggplot(pairwise_data, aes(x=paste(File1, "vs", File2), y=Overlap_Percentage)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x="File Pairs", y="Overlap Percentage (%)",
         title="Pairwise Overlap Percentages")
print(p1)

# Three-way overlap visualization if data exists
if (!is.null(three_way_data)) {
    p2 <- ggplot(three_way_data, 
                 aes(x=paste(File1, "vs", File2, "vs", File3), 
                     y=Overlap_Percentage)) +
        geom_bar(stat="identity", fill="darkred") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x="File Trios", y="Overlap Percentage (%)",
             title="Three-way Overlap Percentages")
    print(p2)
}

dev.off()
EOF

# Run R script if R is available
if command -v Rscript &> /dev/null; then
    echo "Creating visualizations..."
    cd "$OUTPUT_DIR/summary" && Rscript plot_overlaps.R
    
    # Install required R packages if missing
    Rscript -e '
    if (!require("VennDiagram")) {
        install.packages("VennDiagram", repos="https://cloud.r-project.org")
    }
    if (!require("gridExtra")) {
        install.packages("gridExtra", repos="https://cloud.r-project.org")
    }
    '
fi

echo "Analysis complete!"
echo "Results written to: $OUTPUT_DIR"
echo "See summary file: $SUMMARY_FILE"
echo "See pairwise CSV summary: $PAIRWISE_CSV"
if [ $NUM_FILES -ge 3 ]; then
    echo "See three-way CSV summary: $THREE_WAY_CSV"
fi
if [ -f "$OUTPUT_DIR/summary/overlap_visualization.pdf" ]; then
    echo "See visualizations: $OUTPUT_DIR/summary/overlap_visualization.pdf"
fi 