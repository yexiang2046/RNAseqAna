#!/usr/bin/env bash

# Help message
usage() {
    echo "Usage: $0 [-h] -b BED_FILE -f FASTA_FILE -o OUTPUT_FILE"
    echo "Extract RNA sequences from a reference fasta file using bed coordinates"
    echo ""
    echo "Arguments:"
    echo "  -b BED_FILE     Input BED file with coordinates (must be BED6 format with strand in column 6)"
    echo "  -f FASTA_FILE   Reference genome fasta file"
    echo "  -o OUTPUT_FILE  Output fasta file"
    echo "  -h             Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "hb:f:o:" opt; do
    case $opt in
        h) usage ;;
        b) BED_FILE="$OPTARG" ;;
        f) FASTA_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        ?) usage ;;
    esac
done

# Check required arguments
if [ -z "$BED_FILE" ] || [ -z "$FASTA_FILE" ] || [ -z "$OUTPUT_FILE" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input files exist
if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file does not exist: $BED_FILE"
    exit 1
fi

if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: FASTA file does not exist: $FASTA_FILE"
    exit 1
fi

# Check if bedtools is installed
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is not installed"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

echo "Extracting RNA sequences..."
echo "Input BED file: $BED_FILE"
echo "Reference FASTA: $FASTA_FILE"
echo "Output file: $OUTPUT_FILE"

# Check if fasta index exists, if not create it
if [ ! -f "${FASTA_FILE}.fai" ]; then
    echo "Creating FASTA index..."
    samtools faidx "$FASTA_FILE"
fi

# Extract sequences using bedtools getfasta
# -s: Force strandedness (reverse complement if strand is -)
# -name: Use the name field for the fasta header
bedtools getfasta -fi "$FASTA_FILE" \
                 -bed "$BED_FILE" \
                 -fo "$OUTPUT_FILE" \
                 -s \
                 -name

# replace the T with U in the sequence in the output file, without changing the header
sed -i 's/T/U/g' "$OUTPUT_FILE"


# Check if output was created successfully
if [ ! -s "$OUTPUT_FILE" ]; then
    echo "Error: Failed to create output file or output file is empty"
    exit 1
fi

# Print summary
echo "Done!"
echo "Number of sequences extracted: $(grep -c "^>" "$OUTPUT_FILE")"
echo "First few sequences:"
head -n 6 "$OUTPUT_FILE" 