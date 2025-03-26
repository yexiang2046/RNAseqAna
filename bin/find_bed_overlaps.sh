#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i BED_DIR -o OUTPUT_DIR [-m MIN_OVERLAP] [-s STRAND_AWARE] [-H]"
    echo "Find overlapping regions between BED files"
    echo ""
    echo "Arguments:"
    echo "  -i BED_DIR      Directory containing BED files"
    echo "  -o OUTPUT_DIR   Output directory"
    echo "  -m MIN_OVERLAP  Minimum overlap required (default: 1bp)"
    echo "  -s STRAND_AWARE Consider strand in bedtools intersect: true/false (default: true)"
    echo "  -H             Use only host regions (chromosomes starting with 'chr')"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hi:o:m:s:H" opt; do
    case $opt in
        h) usage ;;
        i) BED_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) MIN_OVERLAP="$OPTARG" ;;
        s) STRAND_AWARE="$OPTARG" ;;
        H) HOST_ONLY=true ;;
        ?) usage ;;
    esac
done

# Check required arguments and set defaults
if [ -z "$BED_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

MIN_OVERLAP=${MIN_OVERLAP:-1}
STRAND_AWARE=${STRAND_AWARE:-"true"}
HOST_ONLY=${HOST_ONLY:-false}

# Validate inputs
if [ ! -d "$BED_DIR" ]; then
    echo "Error: BED directory does not exist: $BED_DIR"
    exit 1
fi

if [ ! -z "$STRAND_AWARE" ] && [ "$STRAND_AWARE" != "true" ] && [ "$STRAND_AWARE" != "false" ]; then
    echo "Error: Invalid strand_aware value. Must be 'true' or 'false'"
    exit 1
fi

if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR/pairwise" "$OUTPUT_DIR/three_way" "$OUTPUT_DIR/summary" "$OUTPUT_DIR/filtered"

# Helper functions
get_basename() {
    basename "$1" .bed
}

count_regions() {
    local file="$1"
    local filter="$2"
    if [ "$filter" = "host" ]; then
        awk '$1 ~ /^chr/' "$file" | wc -l
    else
        wc -l < "$file"
    fi
}

filter_host_regions() {
    local input="$1"
    local output="$2"
    awk '$1 ~ /^chr/' "$input" > "$output"
}

# Get list of BED files and filter if needed
if [ "$HOST_ONLY" = true ]; then
    echo "Host-only mode enabled: filtering for chromosomes starting with 'chr'"
    for file in "$BED_DIR"/*.bed; do
        filtered="$OUTPUT_DIR/filtered/$(basename "$file")"
        filter_host_regions "$file" "$filtered"
        BED_FILES+=("$filtered")
    done
else
    BED_FILES=("$BED_DIR"/*.bed)
fi

NUM_FILES=${#BED_FILES[@]}
if [ $NUM_FILES -lt 2 ]; then
    echo "Error: At least 2 BED files are required. Found: $NUM_FILES"
    exit 1
fi

echo "Found $NUM_FILES BED files"
for file in "${BED_FILES[@]}"; do
    echo "File: $file"
    head -n 1 "$file"
done

# Initialize output files
OVERLAP_CSV="$OUTPUT_DIR/summary/overlap_summary.csv"
SUMMARY_FILE="$OUTPUT_DIR/summary/overlap_summary.txt"

# Initialize variables for the main trio
total1=0
total2=0
total3=0
overlap12=0
overlap13=0
overlap23=0
overlap123=0

# Write headers
echo "total1,total2,total3,overlap12,overlap13,overlap23,overlap123" > "$OVERLAP_CSV"
{
    echo "Overlap Analysis Summary"
    echo "Date: $(date)"
    echo "Minimum overlap required: ${MIN_OVERLAP}bp"
    echo "Strand-aware: $STRAND_AWARE"
    echo "Host-only mode: $HOST_ONLY"
    echo "----------------------------------------"
} > "$SUMMARY_FILE"

# Get counts for first three files
echo "Getting counts for first three files..."
for ((i=0; i<3 && i<$NUM_FILES; i++)); do
    file="${BED_FILES[$i]}"
    count=$(count_regions "$file" "$([ "$HOST_ONLY" = true ] && echo "host" || echo "all")")
    echo "File $((i+1)): $count regions"
    case $i in
        0) total1=$count ;;
        1) total2=$count ;;
        2) total3=$count ;;
    esac
done

# Process pairwise overlaps
echo "Processing pairwise overlaps..."
for ((i=0; i<$NUM_FILES; i++)); do
    for ((j=i+1; j<$NUM_FILES; j++)); do
        file1="${BED_FILES[$i]}"
        file2="${BED_FILES[$j]}"
        base1=$(get_basename "$file1")
        base2=$(get_basename "$file2")
        
        echo "Processing pair: $base1 vs $base2"
        echo "File1: $file1"
        echo "File2: $file2"
        
        # Find overlaps
        overlap_file="$OUTPUT_DIR/pairwise/${base1}_vs_${base2}_overlap.bed"
        echo "Running bedtools intersect..."
        if [ "$STRAND_AWARE" = "true" ]; then
            bedtools intersect -a "$file1" -b "$file2" -wa -f "$MIN_OVERLAP" -s | sort -k1,1 -k2,2n -k3,3n | uniq > "$overlap_file"
        else
            bedtools intersect -a "$file1" -b "$file2" -wa -f "$MIN_OVERLAP" | sort -k1,1 -k2,2n -k3,3n | uniq > "$overlap_file"
        fi
        
        # Count overlaps
        overlaps=$(wc -l < "$overlap_file")
        echo "Found $overlaps overlapping regions"
        
        # Store overlaps for main trio
        if [ $i -eq 0 ] && [ $j -eq 1 ]; then
            overlap12=$overlaps
        elif [ $i -eq 0 ] && [ $j -eq 2 ]; then
            overlap13=$overlaps
        elif [ $i -eq 1 ] && [ $j -eq 2 ]; then
            overlap23=$overlaps
        fi
        
        # Write to summary
        {
            echo -e "\nPairwise: $base1 vs $base2:"
            echo "File 1 ($base1): $(count_regions "$file1" "all") regions"
            echo "File 2 ($base2): $(count_regions "$file2" "all") regions"
            echo "Overlaps: $overlaps regions"
        } >> "$SUMMARY_FILE"
    done
done

# Process three-way overlaps if we have at least 3 files
if [ $NUM_FILES -ge 3 ]; then
    echo "Processing three-way overlaps..."
    file1="${BED_FILES[0]}"
    file2="${BED_FILES[1]}"
    file3="${BED_FILES[2]}"
    base1=$(get_basename "$file1")
    base2=$(get_basename "$file2")
    base3=$(get_basename "$file3")
    
    echo "Processing trio: $base1 vs $base2 vs $base3"
    
    # Find three-way overlaps
    temp_overlap="$OUTPUT_DIR/three_way/temp_overlap.bed"
    overlap_file="$OUTPUT_DIR/three_way/${base1}_${base2}_${base3}_overlap.bed"
    
    # First find overlap between first two files
    echo "Running first bedtools intersect..."
    if [ "$STRAND_AWARE" = "true" ]; then
        bedtools intersect -a "$file1" -b "$file2" -wa -f "$MIN_OVERLAP" -s | sort -k1,1 -k2,2n -k3,3n | uniq > "$temp_overlap"
    else
        bedtools intersect -a "$file1" -b "$file2" -wa -f "$MIN_OVERLAP" | sort -k1,1 -k2,2n -k3,3n | uniq > "$temp_overlap"
    fi
    
    # Then find overlap with third file
    echo "Running second bedtools intersect..."
    if [ "$STRAND_AWARE" = "true" ]; then
        bedtools intersect -a "$temp_overlap" -b "$file3" -wa -f "$MIN_OVERLAP" -s | sort -k1,1 -k2,2n -k3,3n | uniq > "$overlap_file"
    else
        bedtools intersect -a "$temp_overlap" -b "$file3" -wa -f "$MIN_OVERLAP" | sort -k1,1 -k2,2n -k3,3n | uniq > "$overlap_file"
    fi
    
    rm -f "$temp_overlap"
    
    # Count three-way overlaps
    overlap123=$(wc -l < "$overlap_file")
    echo "Found $overlap123 three-way overlapping regions"
    
    # Write to summary
    {
        echo -e "\nThree-way: $base1 vs $base2 vs $base3:"
        echo "File 1 ($base1): $(count_regions "$file1" "all") regions"
        echo "File 2 ($base2): $(count_regions "$file2" "all") regions"
        echo "File 3 ($base3): $(count_regions "$file3" "all") regions"
        echo "Common regions: $overlap123"
    } >> "$SUMMARY_FILE"
fi

# Write final CSV line with debug output
echo "Debug values:"
echo "total1: $total1"
echo "total2: $total2"
echo "total3: $total3"
echo "overlap12: $overlap12"
echo "overlap13: $overlap13"
echo "overlap23: $overlap23"
echo "overlap123: $overlap123"

echo "$total1,$total2,$total3,$overlap12,$overlap13,$overlap23,$overlap123" >> "$OVERLAP_CSV"

# Generate Venn diagrams
echo "Generating Venn diagrams..."
Rscript bin/plot_overlaps.R "$OVERLAP_CSV"

echo "Analysis complete. Results written to $OUTPUT_DIR/summary/"

