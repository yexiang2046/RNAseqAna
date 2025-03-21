#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -i INPUT_GTF -o OUTPUT_DIR"
    echo "Extract genomic features from GTF file using gtftools"
    echo ""
    echo "Arguments:"
    echo "  -i INPUT_GTF     Input GTF file"
    echo "  -o OUTPUT_DIR    Output directory for BED files"
    echo "  -h              Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hi:o:" opt; do
    case $opt in
        h) usage ;;
        i) INPUT="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$INPUT" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input file exists
if [ ! -f "$INPUT" ]; then
    echo "Error: Input file does not exist: $INPUT"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Change to output directory
cd "$OUTPUT_DIR"

# Extract features using gtftools
echo "Extracting genes..."
gtftools -g genes.bed "$INPUT"
# format to bed6, add chr prefix
awk 'BEGIN{OFS="\t"} {
    print "chr"$1, $2, $3, $5, $6, $4, $7
}' genes.bed > genes.bed

echo "Extracting exons..."
gtftools -m exons.bed "$INPUT"
# format to bed6, add chr prefix
awk 'BEGIN{OFS="\t"} {
    print "chr"$1, $2, $3, $5, $6, $4, $7
}' exons.bed > exons.bed

echo "Extracting introns..."
gtftools -d introns.bed "$INPUT"
# format to bed6, add chr prefix
awk 'BEGIN{OFS="\t"} {
    print "chr"$1, $2, $3, $5, $4, '.'
}' introns.bed > introns.bed

echo "Extracting CDS..."
gtftools -o cds.bed "$INPUT"
# format to bed6, add chr prefix
awk 'BEGIN{OFS="\t"} {
    print "chr"$1, $2, $3, $5, $6, $4, $7
}' cds.bed > cds.bed


echo "Extracting UTRs..."
gtftools -u utrs.bed "$INPUT"
# format to bed6, add chr prefix
awk 'BEGIN{OFS="\t"} {
    print "chr"$1, $2, $3, $5, $6, $4, $7
}' utrs.bed > utrs.bed

# Split UTRs into 5' and 3'
echo "Splitting UTRs..."
awk '$5=="5UTR" {print > "five_prime_utr.bed"} $5=="3UTR" {print > "three_prime_utr.bed"}' utrs.bed
rm utrs.bed



echo "Feature extraction complete!" 