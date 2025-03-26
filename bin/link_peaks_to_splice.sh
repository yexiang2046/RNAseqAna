#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i PEAKS.bed -s SPLICE_Q.tsv -o OUTPUT_DIR [-d DEBUG_DIR] [-x X_MIN] [-X X_MAX]"
    echo "Link SPLICE-Q introns to overlapping peaks"
    echo ""
    echo "Required arguments:"
    echo "  -i PEAKS.bed     Input peaks file in BED format"
    echo "  -s SPLICE_Q.tsv  SPLICE-Q results file"
    echo "  -o OUTPUT_DIR    Output directory for all files"
    echo "Optional arguments:"
    echo "  -d DEBUG_DIR     Directory to save debug files (default: debug)"
    echo "  -x X_MIN        Minimum value for ridge plot x-axis (default: none)"
    echo "  -X X_MAX        Maximum value for ridge plot x-axis (default: none)"
    echo "  -h              Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hi:s:o:d:x:X:" opt; do
    case $opt in
        h) usage ;;
        i) PEAKS="$OPTARG" ;;
        s) SPLICE_Q="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) DEBUG_DIR="$OPTARG" ;;
        x) X_MIN="$OPTARG" ;;
        X) X_MAX="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Set default debug directory
DEBUG_DIR=${DEBUG_DIR:-"debug"}

# Check required arguments
if [ -z "$PEAKS" ] || [ -z "$SPLICE_Q" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input files exist
if [ ! -f "$PEAKS" ]; then
    echo "Error: Peaks file not found: $PEAKS"
    exit 1
fi

if [ ! -f "$SPLICE_Q" ]; then
    echo "Error: SPLICE-Q file not found: $SPLICE_Q"
    exit 1
fi

# Create output and debug directories
if [ -f "$OUTPUT_DIR" ]; then
    echo "Error: Output directory '$OUTPUT_DIR' exists as a file"
    exit 1
fi
mkdir -p "$OUTPUT_DIR"
mkdir -p "$DEBUG_DIR"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

# Convert SPLICE-Q output to BED format for introns
# Note: For introns, use 5' end as start and 3' end as end
# Coordinates are 1-based in SPLICE-Q output, convert to 0-based
awk -F'\t' 'NR>1 {
    chr=$1
    if (!match(chr, /^chr/)) chr="chr"chr
    # For + strand, 5sj end to 3sj start
    if ($2 == "+") {
        start=$7-1  # Convert to 0-based
        end=$10-1   # Convert to 0-based
    }
    # For - strand, 3sj end to 5sj start
    else {
        start=$11-1  # Convert to 0-based
        end=$6-1     # Convert to 0-based
    }
    # Skip if start >= end (invalid intron)
    if (start >= end) next
    name=$3"_"$4"_"$5
    strand=$2
    score=$14
    # Print original SPLICE-Q fields for debugging
    print $0 > "'$DEBUG_DIR'/splice_q_filtered.tsv"
    # Print transformed intron coordinates (BED6 + score)
    print chr"\t"start"\t"end"\t"name"\t"0"\t"strand"\t"score
}' "$SPLICE_Q" > "$TMP_DIR/introns.bed"

# Copy introns BED file to debug directory
cp "$TMP_DIR/introns.bed" "$DEBUG_DIR/introns.unsorted.bed"

# Sort introns
sort -k1,1 -k2,2n "$TMP_DIR/introns.bed" > "$TMP_DIR/introns.sorted.bed"

# Copy sorted introns to debug directory
cp "$TMP_DIR/introns.sorted.bed" "$DEBUG_DIR/introns.sorted.bed"

# Intersect introns with peaks and extract desired columns
bedtools intersect -a "$TMP_DIR/introns.sorted.bed" -b "$PEAKS" -wa -wb | \
    tee "$DEBUG_DIR/introns_peaks_intersect.bed" | \
    awk -F'\t' 'BEGIN {OFS="\t"} {
        # Keep intron info (chr, start, end, name, score, strand, spliceq_score)
        # and peak info (score, strand)
        print $1, $2, $3, $4, $5, $6, $7, $11, $12
    }' | sort -k1,1 -k2,2n > "$TMP_DIR/introns_with_peaks.bed"

# Export introns without peak overlap to a separate file
bedtools intersect -a "$TMP_DIR/introns.sorted.bed" -b "$PEAKS" -v | \
    tee "$DEBUG_DIR/introns_without_peaks.bed" | \
    sort -k1,1 -k2,2n > "$OUTPUT_DIR/introns_no_overlap.bed"

# Sort final output
sort -k1,1 -k2,2n "$TMP_DIR/introns_with_peaks.bed" > "$OUTPUT_DIR/introns_with_peaks.bed"

# Copy final output to debug directory
cp "$OUTPUT_DIR/introns_with_peaks.bed" "$DEBUG_DIR/final_annotated_introns.bed"

# Generate violin plot of SPLICE-Q scores
echo "Generating violin plot of SPLICE-Q scores..."
if [ ! -f "$OUTPUT_DIR/introns_with_peaks.bed" ] || [ ! -f "$OUTPUT_DIR/introns_no_overlap.bed" ]; then
    echo "Error: Required input files for violin plot not found"
    exit 1
fi

# Build R script command with optional x-axis range parameters
R_CMD="Rscript bin/plot_spliceq_scores.R \"$OUTPUT_DIR/introns_with_peaks.bed\" \"$OUTPUT_DIR/introns_no_overlap.bed\""
if [ ! -z "$X_MIN" ] && [ ! -z "$X_MAX" ]; then
    R_CMD="$R_CMD $X_MIN $X_MAX"
fi

# Run R script
eval "$R_CMD"

# Move the plot files to output directory
if [ -f "spliceq_scores_violin.pdf" ]; then
    mv spliceq_scores_violin.pdf "$OUTPUT_DIR/spliceq_scores_violin.pdf"
else
    echo "Warning: Violin plot file not generated"
fi

if [ -f "spliceq_scores_ridge.pdf" ]; then
    mv spliceq_scores_ridge.pdf "$OUTPUT_DIR/spliceq_scores_ridge.pdf"
else
    echo "Warning: Ridge plot file not generated"
fi

echo "Debug files written to: $DEBUG_DIR/"
echo "Files:"
echo "  - splice_q_filtered.tsv: Original SPLICE-Q entries used"
echo "  - introns.unsorted.bed: Transformed intron coordinates (unsorted)"
echo "  - introns.sorted.bed: Transformed intron coordinates (sorted, BED6+score format)"
echo "  - introns_peaks_intersect.bed: Raw bedtools intersect output"
echo "  - introns_without_peaks.bed: Introns with no peak overlap"
echo "  - final_annotated_introns.bed: Final output with peak info"
echo ""
echo "Output files in $OUTPUT_DIR/:"
echo "  - introns_with_peaks.bed: Introns with peak overlap"
echo "  - introns_no_overlap.bed: Introns without peak overlap"
echo "  - spliceq_scores_violin.pdf: Violin plot of SPLICE-Q scores"
echo "  - spliceq_scores_ridge.pdf: Ridge plot of SPLICE-Q scores"
echo ""
echo "Output format: chr start end intron_id score strand spliceq_score peak_score peak_strand"
echo "Done! All files written to: $OUTPUT_DIR/" 