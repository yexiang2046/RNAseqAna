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

echo "Extracting exons..."
gtftools -m exons.bed "$INPUT"

echo "Extracting introns..."
gtftools -d introns.bed "$INPUT"

echo "Extracting CDS..."
gtftools -o cds.bed "$INPUT"

echo "Extracting splice regions..."
gtftools -q splice_regions.bed "$INPUT"

echo "Extracting UTRs..."
gtftools -u utrs.bed "$INPUT"

# Split UTRs into 5' and 3'
echo "Splitting UTRs..."
awk '$6=="5UTR" {print > "five_prime_utr.bed"} $6=="3UTR" {print > "three_prime_utr.bed"}' utrs.bed
rm utrs.bed

# Add headers to BED files
echo "Adding headers to BED files..."
for bed in *.bed; do
    echo -e "chr\tstart\tend\tname\tscore\tstrand\tgene_name\tgene_id" > "tmp_${bed}"
    cat "$bed" >> "tmp_${bed}"
    mv "tmp_${bed}" "$bed"
done

echo "Feature extraction complete!" 